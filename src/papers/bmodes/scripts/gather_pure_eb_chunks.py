"""Gather MC sample chunks and compute the final pure E/B covariance.

CLI refactor of the former Snakemake ``script:`` gather rule. Reads the actual
ξ± data vectors (reporting + integration grids), computes the pure E/B/ambiguous
decomposition (Schneider 2022) in the geometry ``--integration-method``
implies — pointwise on the reporting grid under ``adaptive_dqags``, on the fine
grid then averaged into the reporting bins under ``fixed_operator``, with the
weights ``--rebinning-weights`` selects — stacks the per-chunk MC sample blocks, and
forms the 6-block covariance under ``--covariance-method``: the empirical
covariance of the Monte Carlo chunks, or the exact ``M C_ξ Mᵀ`` of the
averaged transform (no chunks; ``--cov-integration`` is C_ξ). Writes the
per-(version, blind) ``<version>_<blind>_pure_eb_semianalytic.npz`` consumed
by every downstream pure-mode plot / PTE, recording which construction made
the covariance (``pure_mode_covariance.covariance_record``).

    python gather_pure_eb_chunks.py \
        --version SP_v1.4.6.3_leak_corr --blind A \
        --xi-reporting  <xi 20-bin .txt> \
        --xi-integration <xi 1000-bin .txt> \
        --chunks-dir <dir with pure_eb_chunk_*.npz> \
        --min-sep 1.0 --max-sep 250.0 --nbins 20 \
        --min-sep-int 0.5 --max-sep-int 300.0 --nbins-int 1000 \
        --integration-method adaptive_dqags --rebinning-weights none \
        --covariance-method monte_carlo --out <output_dir>
"""

import argparse
import glob
import os

import numpy as np
from pure_eb_transform import (
    AveragedPureEBTransform,
    PureEBTransform,
    add_integration_limits_argument,
    add_integration_method_argument,
    add_rebinning_weights_argument,
    check_rebinning_weights,
    geometry_for,
    resolve_integration_limits,
)
from pure_mode_covariance import (
    add_covariance_method_argument,
    covariance_record,
    exact_covariance,
)


def _load_xi(path, nbins):
    """Load ξ± from a TreeCorr text dump."""
    data = np.loadtxt(path, comments="#", max_rows=nbins)
    loaded = {"meanr": data[:, 1], "xip": data[:, 3], "xim": data[:, 4]}
    # TreeCorr's ASCII header is
    # r_nom meanr meanlogr xip xim xip_im xim_im sigma_xip sigma_xim weight npairs
    if data.shape[1] > 10:
        loaded["npairs"] = data[:, 10]
    return loaded


def data_vector(
    gg,
    gg_int,
    *,
    min_sep,
    max_sep,
    nbins,
    integration_method,
    rebinning_weights,
    integration_limits,
    geometry=None,
    local_from_eval_grid=False,
):
    """The pure-mode DATA vector, from the measured reporting/integration xi.

    Factored out of :func:`gather` so a second caller can build the same
    vector without also gathering MC chunks: ``pure_eb_scan_statistics`` needs
    it for a variant whose covariance is aliased to another's, since the
    borrowed bundle carries the donor's data vector rather than that variant's.
    Splitting it rather than reimplementing it is what makes the two paths the
    same by construction; the caller checks that against a stored vector.

    *gg* / *gg_int* are the dicts :func:`_load_xi` returns.  Returns the
    six-tuple in ``REFERENCE_ORDER``.
    """

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
        local_from_eval_grid=local_from_eval_grid,
    )
    if isinstance(transform, AveragedPureEBTransform):
        return transform(gg_int["xip"], gg_int["xim"], gg["xip"], gg["xim"])
    return transform(gg["xip"], gg["xim"], gg_int["xip"], gg_int["xim"])


def build_transform(
    gg,
    gg_int,
    *,
    min_sep,
    max_sep,
    nbins,
    integration_method,
    rebinning_weights,
    integration_limits,
    geometry=None,
    local_from_eval_grid=False,
):
    """The pure-mode transform the decisions select, built on the data grids."""
    geometry = geometry_for(integration_method, geometry)
    theta_window = resolve_integration_limits(integration_limits, geometry)
    check_rebinning_weights(rebinning_weights, geometry)
    if geometry == "averaged":
        if rebinning_weights == "npairs" and "npairs" not in gg_int:
            raise ValueError(
                "the integration xi has no npairs column; --rebinning-weights "
                "npairs weights the fine-to-reporting average by pair counts"
            )
        edges = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins + 1)
        return AveragedPureEBTransform(
            integration_method,
            left_edges=edges[:-1],
            right_edges=edges[1:],
            theta_int=gg_int["meanr"],
            rebinning_weights=rebinning_weights,
            npairs=gg_int.get("npairs"),
            theta_window=theta_window,
            local_from_eval_grid=local_from_eval_grid,
        )
    return PureEBTransform(
        integration_method,
        theta=gg["meanr"],
        theta_int=gg_int["meanr"],
        tmin=min_sep,
        tmax=max_sep,
        # Mirror precompute_pure_eb_chunk: the single-grid operator evaluates
        # at nodes, so the pointwise operator chain snaps to the nearest node.
        snap_eval=integration_method == "fixed_operator",
        # The hybrid corner: the integral term is still the snapped
        # operator's, but the local 0.5*(xi_p +/- xi_m) term is read from
        # the reporting-grid draw rather than from one fine node.
        local_from_eval_grid=local_from_eval_grid,
    )


def gather(
    version,
    blind,
    xi_reporting,
    xi_integration,
    chunk_files,
    min_sep,
    max_sep,
    nbins,
    nbins_int,
    output_dir,
    integration_method,
    rebinning_weights,
    integration_limits,
    covariance_method,
    cov_integration=None,
    geometry=None,
    local_from_eval_grid=False,
):
    print(f"Gathering pure E/B for blind {blind} ({covariance_method} covariance)")

    gg = _load_xi(xi_reporting, nbins)
    gg_int = _load_xi(xi_integration, nbins_int)

    # geometry=None keeps the method's default mapping (see precompute_pure_eb_chunk).
    eb_results = data_vector(
        gg,
        gg_int,
        min_sep=min_sep,
        max_sep=max_sep,
        nbins=nbins,
        integration_method=integration_method,
        rebinning_weights=rebinning_weights,
        integration_limits=integration_limits,
        geometry=geometry,
        local_from_eval_grid=local_from_eval_grid,
    )
    xip_E, xim_E, xip_B, xim_B, xip_amb, xim_amb = eb_results

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
        if local_from_eval_grid or not isinstance(transform, AveragedPureEBTransform):
            raise ValueError(
                "exact_linear needs the averaged transform on the fine grid alone "
                "(fixed_operator, no --local-from-eval-grid)"
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

    package = {
        "theta": gg["meanr"],
        "theta_int": gg_int["meanr"],
        "xip_total": gg["xip"],
        "xim_total": gg["xim"],
        "xip_E": xip_E,
        "xim_E": xim_E,
        "xip_B": xip_B,
        "xim_B": xim_B,
        "xip_amb": xip_amb,
        "xim_amb": xim_amb,
        "cov_pure_eb": cov_pure_eb,
        **record,
    }

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{version}_{blind}_pure_eb_semianalytic.npz")
    np.savez(out_path, **package)
    print(f"Saved to {out_path}")
    return out_path


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
    ap.add_argument("--version", required=True)
    ap.add_argument("--blind", default="A")
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
    ap.add_argument("--npatch", type=int, default=1)
    ap.add_argument("--out", required=True, help="Output directory (lc {output})")
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
    ap.add_argument(
        "--local-from-eval-grid",
        action="store_true",
        help='Pointwise fixed_operator only: take the local 0.5*(xi_p +/- xi_m) term from the reporting-grid xi at the requested abscissa rather than from the snapped fine node, as the dqags pointwise chain does. The integral term still comes from the fine grid at the snapped node.',
    )
    a = ap.parse_args(argv)
    exact = a.covariance_method == "exact_linear"
    if exact and not a.cov_integration:
        ap.error("--covariance-method exact_linear needs --cov-integration")
    gather(
        version=a.version,
        blind=a.blind,
        xi_reporting=a.xi_reporting,
        xi_integration=a.xi_integration,
        chunk_files=None if exact else _resolve_chunks(a),
        min_sep=a.min_sep,
        max_sep=a.max_sep,
        nbins=a.nbins,
        nbins_int=a.nbins_int,
        output_dir=a.out,
        integration_method=a.integration_method,
        rebinning_weights=a.rebinning_weights,
        integration_limits=a.integration_limits,
        covariance_method=a.covariance_method,
        cov_integration=a.cov_integration,
        geometry=a.geometry,
        local_from_eval_grid=a.local_from_eval_grid,
    )


if __name__ == "__main__":
    _from_cli()
