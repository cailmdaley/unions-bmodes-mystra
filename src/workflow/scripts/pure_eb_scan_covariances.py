"""Both pure-mode covariances for ONE decision variant, in one bundle.

The scan's only heavy stage, and one output per variant so the four run in
parallel on separate nodes rather than serially inside a single job.  It
produces the two MC-propagated pure-mode covariances the comparison needs:

* ``<variant>__data.npz`` -- the data's own pure-mode vectors, from
  ``pure_eb_modes.py``, with their 120x120 covariance, from
  ``run_pure_eb_covariance.sh``.  The artifact the variant's DATA PTE is
  computed from downstream.
* ``<variant>__mockmatched.npz`` -- the same propagation at the mock suite's
  own noise, from ``mock_matched_pure_eb_cov.py``.  The twin of
  ``pure_cov_mockmatched``.

Both are needed *per variant*, not shared.  The transform that propagates the
2000 draws is the thing under scan: dqags' quadrature error is noise per draw
and inflates the propagated variance (xi_-^B bin 0 by x3.3, xi_-^amb by x17),
so reusing one covariance across variants would compare each variant's data
vector against another variant's error bars.  The seeds are identical across
variants (42 + chunk/batch id, same block sizes), so the variant-to-variant
difference carries no Monte-Carlo noise of its own.

The scan runs 500 draws, not the 2000 the paper's fiducial covariance uses.
At 2000 the four dqags-averaged propagations alone would have run for days,
and the scan does not need the fiducial's per-element precision: it compares
variants against each other, and 500 draws give ~0.4% per-element noise
against the 10--30% effects being resolved.  The fiducial 2000-draw
covariance is a separate artifact and is untouched by this stage.

    python pure_eb_scan_covariances.py --xi-reporting R --xi-integration I \
        --cov-integration C --cov-integration-mockmatched M \
        --variant averaged_uniform_fixedop --out-file BUNDLE.tar.gz
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np

from pure_eb_scan_variants import BY_NAME, NAMES

# Mirrors pure_eb_transform.GEOMETRIES. A variant running in its method's
# default geometry passes no --geometry at all, so its two chains stay
# byte-for-byte the calls the single-variant outputs make.
DEFAULT_GEOMETRY = {"adaptive_dqags": "pointwise", "fixed_operator": "averaged"}

FIDUCIAL_VERSION = "SP_v1.4.6.3_leak_corr"
FIDUCIAL_BLIND = "A"
CAT_CONFIG = "src/cosmo_val/cat_config.yaml"


def _run(command: list[str], tag: str) -> None:
    print(f"\n=== {tag} ===\n{' '.join(command)}", flush=True)
    started = time.perf_counter()
    result = subprocess.run(command)
    if result.returncode != 0:
        raise SystemExit(
            f"{tag} failed with exit code {result.returncode}; the scan stops "
            "rather than bundling a partial set of covariances."
        )
    print(f"=== {tag}: {time.perf_counter() - started:.1f}s ===", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xi-reporting", required=True)
    parser.add_argument("--xi-integration", required=True)
    parser.add_argument("--cov-integration", required=True)
    parser.add_argument("--cov-integration-mockmatched", required=True)
    parser.add_argument(
        "--n-samples",
        type=int,
        default=500,
        help=(
            "MC draws propagated through the transform. The scan's 500 is not "
            "the fiducial covariance's 2000: see the module docstring."
        ),
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-file", required=True)
    parser.add_argument(
        "--variant",
        required=True,
        choices=list(NAMES),
        help="The single scan variant this invocation builds covariances for.",
    )
    args = parser.parse_args()

    variants = [BY_NAME[args.variant]]
    aliased = [v for v in variants if v.covariance_from]
    if aliased:
        variant = aliased[0]
        raise SystemExit(
            f"{variant.name} declares covariance_from="
            f"{variant.covariance_from!r}: its propagated covariance is the "
            f"same object as {variant.covariance_from}'s, so it has no bundle "
            "of its own and no output should ask for one. Materialize "
            f"pure_eb_scan_covariance_{variant.covariance_from} instead; "
            "pure_eb_scan_statistics resolves the alias."
        )
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        staging = Path(tmp)
        provenance: dict[str, object] = {}
        for variant in variants:
            # The geometry reaches the two chains through the flag each one
            # already carries: --geometry on the shell driver (which forwards
            # it to precompute_pure_eb_chunk) and --transform-averaged on the
            # mock-matched twin. Omitted entirely for a variant that runs in
            # its method's default geometry, so those calls are byte-for-byte
            # the ones the single-variant outputs make.
            geometry_flag: list[str] = []
            if variant.geometry != DEFAULT_GEOMETRY[variant.integration_method]:
                geometry_flag = ["--geometry", variant.geometry]
            # The hybrid pointwise operator variant differs from its neighbour
            # only here: same snap, same integral term, local term taken from
            # the reporting draw. Both chains carry the same flag name.
            local_flag: list[str] = (
                ["--local-from-eval-grid"] if variant.local_from_eval_grid else []
            )

            data_out = staging / f"{variant.name}__data.npz"
            modes_out = staging / f"{variant.name}__modes.npz"
            _run(
                [
                    sys.executable,
                    "src/papers/bmodes/scripts/pure_eb_modes.py",
                    "--xi-reporting", args.xi_reporting,
                    "--xi-integration", args.xi_integration,
                    "--integration-method", variant.integration_method,
                    "--integration-limits", variant.integration_limits,
                    "--rebinning-weights", variant.rebinning_weights,
                    *geometry_flag,
                    *local_flag,
                    "--out-file", str(modes_out),
                ],
                f"{variant.name} / data vector",
            )
            _run(
                [
                    "bash",
                    "src/papers/bmodes/scripts/run_pure_eb_covariance.sh",
                    "--version", FIDUCIAL_VERSION,
                    "--blind", FIDUCIAL_BLIND,
                    "--cat-config", CAT_CONFIG,
                    "--xi-reporting", args.xi_reporting,
                    "--xi-integration", args.xi_integration,
                    "--cov-integration", args.cov_integration,
                    "--integration-method", variant.integration_method,
                    "--integration-limits", variant.integration_limits,
                    "--rebinning-weights", variant.rebinning_weights,
                    *geometry_flag,
                    *local_flag,
                    # The driver defaults to the fiducial 2000; the scan is 500.
                    # The scan's covariances are the Monte Carlo propagation it compares.
                    "--covariance-method", "monte_carlo",
                    "--n-samples", str(args.n_samples),
                    "--out-file", str(data_out),
                ],
                f"{variant.name} / data covariance",
            )
            with np.load(modes_out) as modes, np.load(data_out) as covariance:
                bundle = {**modes, **covariance}
            np.savez(data_out, **bundle)
            modes_out.unlink()

            mock_out = staging / f"{variant.name}__mockmatched.npz"
            _run(
                [
                    sys.executable,
                    "src/workflow/scripts/mock_matched_pure_eb_cov.py",
                    "--cov-integration", args.cov_integration_mockmatched,
                    "--xi-reporting", args.xi_reporting,
                    "--xi-integration", args.xi_integration,
                    "--version", FIDUCIAL_VERSION,
                    "--blind", FIDUCIAL_BLIND,
                    "--cat-config", CAT_CONFIG,
                    "--integration-method", variant.integration_method,
                    "--integration-limits", variant.integration_limits,
                    "--rebinning-weights", variant.rebinning_weights,
                    *geometry_flag,
                    *local_flag,
                    # The scan's covariances are the Monte Carlo propagation it compares.
                    "--covariance-method", "monte_carlo",
                    "--n-samples", str(args.n_samples),
                    "--seed", str(args.seed),
                    "--out", str(mock_out),
                ],
                f"{variant.name} / mock-matched covariance",
            )
            provenance[variant.name] = {
                **variant.as_dict(),
                "data": data_out.name,
                "mockmatched": mock_out.name,
            }

        (staging / "provenance.json").write_text(
            json.dumps(
                {
                    "n_samples": args.n_samples,
                    "seed": args.seed,
                    "version": FIDUCIAL_VERSION,
                    "blind": FIDUCIAL_BLIND,
                    "variants": provenance,
                },
                indent=2,
            )
            + "\n"
        )
        with tarfile.open(out_path, "w:gz") as bundle:
            for path in sorted(staging.iterdir()):
                bundle.add(path, arcname=path.name)

    print(
        f"Saved {out_path}: {args.variant}, 2 covariances at "
        f"{args.n_samples} draws",
        flush=True,
    )


if __name__ == "__main__":
    main()
