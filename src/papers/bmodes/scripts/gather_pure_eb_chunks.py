"""Gather the pure E/B covariance from MC chunks or the exact linear transform.

Forms the 6-block covariance under ``--covariance-method``: the empirical
covariance of the Monte Carlo chunks, or the exact ``M C_ξ Mᵀ`` of the averaged
transform (no chunks; ``--cov-integration`` is C_ξ). The transform is built on
the measured ξ± grids in the geometry ``--integration-method`` implies, with the
weights ``--rebinning-weights`` selects; only the bin centres and pair counts of
those files enter. Writes the npz every pure-mode plot and PTE reads the
covariance from, recording which construction made it
(``pure_mode_covariance.covariance_record``).

    python gather_pure_eb_chunks.py \
        --xi-reporting  <xi 20-bin .txt> \
        --xi-integration <xi 1000-bin .txt> \
        --chunks-dir <dir with pure_eb_chunk_*.npz> \
        --min-sep 1.0 --max-sep 250.0 --nbins 20 \
        --min-sep-int 0.5 --max-sep-int 300.0 --nbins-int 1000 \
        --integration-method adaptive_dqags --rebinning-weights none \
        --covariance-method monte_carlo --out-file <covariance.npz>
"""

import argparse
import glob
import os

import numpy as np
from pure_eb_transform import (
    AveragedPureEBTransform,
    add_integration_limits_argument,
    add_integration_method_argument,
    add_rebinning_weights_argument,
    build_transform,
    load_xi,
)
from pure_mode_covariance import (
    add_covariance_method_argument,
    covariance_record,
    exact_covariance,
)


def gather(
    xi_reporting,
    xi_integration,
    chunk_files,
    min_sep,
    max_sep,
    nbins,
    nbins_int,
    out_file,
    integration_method,
    rebinning_weights,
    integration_limits,
    covariance_method,
    cov_integration=None,
    geometry=None,
):
    print(f"Gathering the pure E/B covariance ({covariance_method})")

    gg = load_xi(xi_reporting, nbins)
    gg_int = load_xi(xi_integration, nbins_int)

    if covariance_method == "exact_linear":
        transform = build_transform(
            gg,
            gg_int,
            min_sep=min_sep,
            max_sep=max_sep,
            nbins=nbins,
            integration_method=integration_method,
            rebinning_weights=rebinning_weights,
            integration_limits=integration_limits,
            geometry=geometry,
        )
        if not isinstance(transform, AveragedPureEBTransform):
            raise ValueError(
                "exact_linear needs the averaged transform on the fine grid "
                "(fixed_operator)"
            )
        cov_pure_eb = exact_covariance(transform, np.loadtxt(cov_integration))
        record = covariance_record("exact_linear")
    else:
        all_samples = []
        for chunk_file in sorted(chunk_files):
            data = np.load(chunk_file)
            all_samples.append(data["eb_samples"])
            print(f"Loaded {len(data['eb_samples'])} samples from {chunk_file}")

        eb_samples = np.vstack(all_samples)
        print(f"Total samples: {len(eb_samples)}")

        cov_pure_eb = np.cov(eb_samples.T)
        record = covariance_record("monte_carlo", n_samples=len(eb_samples))

    np.savez(
        out_file,
        theta=gg["meanr"],
        theta_int=gg_int["meanr"],
        cov_pure_eb=cov_pure_eb,
        **record,
    )
    print(f"Saved to {out_file}")
    return out_file


def _resolve_chunks(args):
    if args.chunks:
        files = args.chunks
    else:
        files = glob.glob(os.path.join(args.chunks_dir, "pure_eb_chunk_*.npz"))
    if not files:
        raise SystemExit(f"No chunk .npz found (chunks-dir={args.chunks_dir})")
    return files


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xi-reporting", required=True)
    ap.add_argument("--xi-integration", required=True)
    ap.add_argument("--chunks-dir", help="Directory holding pure_eb_chunk_*.npz (monte_carlo)")
    ap.add_argument(
        "--cov-integration", help="Fine ξ_± covariance C_ξ (exact_linear)"
    )
    ap.add_argument("--chunks", nargs="+", help="Explicit chunk .npz paths")
    ap.add_argument("--min-sep", type=float, default=1.0)
    ap.add_argument("--max-sep", type=float, default=250.0)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--min-sep-int", type=float, default=0.5)
    ap.add_argument("--max-sep-int", type=float, default=300.0)
    ap.add_argument("--nbins-int", type=int, default=1000)
    ap.add_argument("--out-file", required=True, help="Output npz (lc {output})")
    add_integration_method_argument(ap)
    add_integration_limits_argument(ap)
    add_rebinning_weights_argument(ap)
    add_covariance_method_argument(ap)
    ap.add_argument(
        "--geometry",
        choices=["pointwise", "averaged"],
        default=None,
        help="Override the geometry the integration method defaults to.",
    )
    a = ap.parse_args(argv)
    exact = a.covariance_method == "exact_linear"
    if exact and not a.cov_integration:
        ap.error("--covariance-method exact_linear needs --cov-integration")
    gather(
        xi_reporting=a.xi_reporting,
        xi_integration=a.xi_integration,
        chunk_files=None if exact else _resolve_chunks(a),
        min_sep=a.min_sep,
        max_sep=a.max_sep,
        nbins=a.nbins,
        nbins_int=a.nbins_int,
        out_file=a.out_file,
        integration_method=a.integration_method,
        rebinning_weights=a.rebinning_weights,
        integration_limits=a.integration_limits,
        covariance_method=a.covariance_method,
        cov_integration=a.cov_integration,
        geometry=a.geometry,
    )


if __name__ == "__main__":
    _from_cli()
