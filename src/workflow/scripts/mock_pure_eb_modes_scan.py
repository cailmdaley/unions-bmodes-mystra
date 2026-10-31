"""Mock pure E/B modes for every decision variant, in one npz.

The scan's counterpart of ``mock_pure_eb_scatter``: the same
``compute_ensemble`` transform, looped over the variants in
``pure_eb_scan_variants`` and written as one npz whose arrays are keyed
``<variant>__<block>``.  One output rather than one per variant, because a
variant is not independently interesting here -- the comparison across them is.

    python mock_pure_eb_modes_scan.py --fine-ensemble F --reporting-ensemble R \
        --variants averaged_uniform_fixedop ... --out-file O
"""

from __future__ import annotations

import argparse
import json
import tempfile
import time
from pathlib import Path

import numpy as np

import mock_pure_eb_scatter as scatter  # importable: see its __main__ guard
from pure_eb_scan_variants import VARIANTS, add_variants_argument, resolve

# The six pure-mode blocks, in the reference order.
BLOCKS = ("xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")
# Carried per variant so a consumer never has to re-derive the geometry.
SCALARS = (
    "tmin",
    "tmax",
    "integration_limits",
    "transform",
    "rebinning_weights",
    # True only for the hybrid pointwise operator variant, whose local term is
    # read from the reporting xi rather than the snapped fine node. Taken from
    # the run rather than the table, so provenance records what actually ran.
    "local_from_eval_grid",
    # Largest |theta_snapped/theta_requested - 1|; nonzero only for the
    # pointwise operator chain, whose evaluation points are snapped onto the
    # integration grid. Recorded so that variant is never read without it.
    "snap_offset_max",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fine-ensemble", required=True)
    parser.add_argument("--reporting-ensemble", required=True)
    parser.add_argument("--out-file", required=True)
    add_variants_argument(parser)
    args = parser.parse_args()

    variants = resolve(args.variants)
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    arrays: dict[str, np.ndarray] = {}
    provenance: dict[str, object] = {}
    shared_written = False

    with tempfile.TemporaryDirectory() as tmp:
        for variant in variants:
            started = time.perf_counter()
            # compute_ensemble writes an npz; it is the same call the
            # single-variant output makes, so the two stay bit-identical.
            per_variant = Path(tmp) / f"{variant.name}.npz"
            scatter.compute_ensemble(
                args.fine_ensemble,
                args.reporting_ensemble,
                per_variant,
                integration_method=variant.integration_method,
                rebinning_weights=variant.rebinning_weights,
                integration_limits=variant.integration_limits,
                geometry=variant.geometry,
                local_from_eval_grid=variant.local_from_eval_grid,
            )
            elapsed = time.perf_counter() - started
            with np.load(per_variant, allow_pickle=False) as data:
                for block in BLOCKS:
                    arrays[f"{variant.name}__{block}"] = np.asarray(data[block])
                arrays[f"{variant.name}__theta"] = np.asarray(data["theta"])
                if not shared_written:
                    # Identical across variants by construction: the reporting
                    # binning and the realisation list do not depend on the
                    # transform. Stored once, and asserted below.
                    for key in ("mock_ids", "reporting_left_edges", "reporting_right_edges"):
                        arrays[key] = np.asarray(data[key])
                    shared_written = True
                else:
                    if not np.array_equal(arrays["mock_ids"], np.asarray(data["mock_ids"])):
                        raise ValueError(
                            f"{variant.name} transformed a different realisation "
                            "set than the first variant"
                        )
                provenance[variant.name] = {
                    **variant.as_dict(),
                    **{key: np.asarray(data[key]).item() for key in SCALARS},
                    "seconds": round(elapsed, 1),
                }
            print(f"{variant.name}: {elapsed:.1f}s", flush=True)

    np.savez(
        out_path,
        variants=np.asarray([variant.name for variant in variants]),
        provenance=np.asarray(json.dumps(provenance)),
        **arrays,
    )
    n_real = arrays["mock_ids"].size
    print(
        f"Saved {out_path}: {len(variants)} variants x {n_real} realisations",
        flush=True,
    )


if __name__ == "__main__":
    main()
