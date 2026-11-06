"""Compute one chunk of MC samples for the pure E/B covariance.

Draw ``n_samples // n_chunks`` Gaussian realisations of ξ±(θ) from
the 1000-bin integration-grid CosmoCov covariance (deterministic seed
``42 + chunk_id``), rebin to the reporting grid, and push each draw through the
Schneider-2022 pure-mode integral transforms (``cosmo_numba``, under the
quadrature ``--integration-method`` selects). ``--integration-method`` also
selects the geometry (``pure_eb_transform.geometry_for``): ``fixed_operator``
skips the input rebinning and instead transforms each draw on the fine grid,
averaging the pure modes into the reporting bins under the weights
``--rebinning-weights`` selects. The per-chunk
E/B/amb sample block is written to ``{out}/pure_eb_chunk_{chunk_id}.npz`` for
the gather stage. Each chunk is independent (fresh RNG per chunk_id), so the
20 chunks give the same 2000-sample covariance whether run in parallel or
looped in one process.

    python precompute_pure_eb_chunk.py \
        --chunk-id 0 --n-chunks 20 --n-samples 2000 \
        --version SP_v1.4.6.3_leak_corr --blind A \
        --cat-config /path/cosmo_val/cat_config.yaml \
        --xi-reporting  <xi 20-bin .txt> \
        --xi-integration <xi 1000-bin .txt> \
        --cov-integration <cov ..._processed.txt> \
        --min-sep 1.0 --max-sep 250.0 --nbins 20 \
        --min-sep-int 0.5 --max-sep-int 300.0 --nbins-int 1000 \
        --npatch 1 --integration-method adaptive_dqags \
        --rebinning-weights none --out <output_dir>
"""

import argparse
import os

import numpy as np
import tqdm
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
from scipy import sparse


def _build_cosmology(cosmo_params):
    """Build a CCL cosmology from a PLANCK18-style params dict."""
    import pyccl as ccl

    return ccl.Cosmology(
        Omega_c=cosmo_params["Omega_m"] - cosmo_params["Omega_b"],
        Omega_b=cosmo_params["Omega_b"],
        h=cosmo_params["h"],
        sigma8=cosmo_params["sigma_8"],
        n_s=cosmo_params["n_s"],
    )


def _load_xi(path, min_sep, max_sep, nbins):
    """Load ξ± from a TreeCorr text dump and recompute the log bin edges."""
    data = np.loadtxt(path, comments="#", max_rows=nbins)
    bin_edges = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins + 1)
    loaded = {
        "meanr": data[:, 1],
        "xip": data[:, 3],
        "xim": data[:, 4],
        "left_edges": bin_edges[:-1],
        "right_edges": bin_edges[1:],
    }
    # TreeCorr's ASCII header is
    # r_nom meanr meanlogr xip xim xip_im xim_im sigma_xip sigma_xim weight npairs
    if data.shape[1] > 10:
        loaded["npairs"] = data[:, 10]
    return loaded


def compute_chunk(
    chunk_id,
    n_chunks,
    n_samples_total,
    version,
    blind,
    cat_config,
    xi_reporting,
    xi_integration,
    cov_integration,
    min_sep,
    max_sep,
    nbins,
    min_sep_int,
    max_sep_int,
    nbins_int,
    output_dir,
    integration_method,
    rebinning_weights,
    integration_limits,
    cosmo_params=None,
    geometry=None,
    local_from_eval_grid=False,
):
    from sp_validation.cosmo_val import CosmologyValidation
    from sp_validation.cosmology import PLANCK18, get_theo_xi

    if cosmo_params is None:
        cosmo_params = dict(PLANCK18)

    # Before the covariance load and the n(z) lookup: a universe that pairs the
    # quadrature and the weighting incoherently should say so immediately.
    # geometry=None keeps the method's default mapping; passing it explicitly
    # is how the decision scan runs adaptive_dqags in the averaged geometry,
    # which is what separates the quadrature effect from the geometry effect.
    geometry = geometry_for(integration_method, geometry)
    theta_window = resolve_integration_limits(integration_limits, geometry)
    check_rebinning_weights(rebinning_weights, geometry)

    samples_per_chunk = n_samples_total // n_chunks
    start_idx = chunk_id * samples_per_chunk
    end_idx = (
        start_idx + samples_per_chunk if chunk_id < n_chunks - 1 else n_samples_total
    )
    n_samples_chunk = end_idx - start_idx
    print(
        f"Chunk {chunk_id}/{n_chunks}: samples {start_idx}-{end_idx} "
        f"({n_samples_chunk} samples)"
    )

    gg = _load_xi(xi_reporting, min_sep, max_sep, nbins)
    gg_int = _load_xi(xi_integration, min_sep_int, max_sep_int, nbins_int)

    cv = CosmologyValidation(
        versions=[version],
        catalog_config=cat_config,
        output_dir=output_dir,
    )
    cv.blind = blind
    z, nz = cv.get_redshift(version)
    z_dist = np.column_stack([z, nz])
    print(f"Using n(z) for blind {blind}")

    cosmo_cov = _build_cosmology(cosmo_params)

    cov_int = np.loadtxt(cov_integration)

    theta_int = gg_int["meanr"]
    # The averaged reporting-local-term variant needs the reporting draws too,
    # so the same uniform fine-to-reporting matrix is built for it.
    if geometry == "pointwise" or local_from_eval_grid:
        reporting_bin_edges = np.concatenate(
            [gg["left_edges"], [gg["right_edges"][-1]]]
        )
        bin_indices = np.digitize(theta_int, reporting_bin_edges) - 1
        valid_mask = (bin_indices >= 0) & (bin_indices < len(gg["meanr"]))
        row_indices, col_indices = bin_indices[valid_mask], np.where(valid_mask)[0]

        binning_matrix = sparse.csr_matrix(
            (np.ones(len(row_indices)), (row_indices, col_indices)),
            shape=(len(gg["meanr"]), nbins_int),
        )
        row_sums = np.array(binning_matrix.sum(axis=1)).flatten()
        binning_matrix = sparse.diags(1 / row_sums) @ binning_matrix

    mean_int = np.concatenate(
        get_theo_xi(
            theta=theta_int,
            z=z_dist[:, 0],
            nz=z_dist[:, 1],
            backend="ccl",
            cosmo=cosmo_cov,
        )
    )

    rng = np.random.default_rng(seed=42 + chunk_id)

    samples_int = rng.multivariate_normal(mean_int, cov_int, size=n_samples_chunk)
    samples_int_xip = samples_int[:, :nbins_int]
    samples_int_xim = samples_int[:, nbins_int:]

    if geometry == "averaged":
        if rebinning_weights == "npairs" and "npairs" not in gg_int:
            raise ValueError(
                f"{xi_integration} has no npairs column; --rebinning-weights "
                "npairs weights the fine-to-reporting average by pair counts"
            )
        transform = AveragedPureEBTransform(
            integration_method,
            left_edges=gg["left_edges"],
            right_edges=gg["right_edges"],
            theta_int=theta_int,
            rebinning_weights=rebinning_weights,
            npairs=gg_int.get("npairs"),
            theta_window=theta_window,
            local_from_eval_grid=local_from_eval_grid,
        )
        if local_from_eval_grid:
            # NOTE these reporting draws are the fine draws uniformly rebinned,
            # so under `uniform` weights they equal the local term the averaged
            # transform would have built anyway: the substitution is a no-op on
            # the MC draws and bites only on the measured data vector, whose
            # reporting xi is its own pair-weighted TreeCorr measurement.
            samples_rep_xip = (binning_matrix @ samples_int_xip.T).T
            samples_rep_xim = (binning_matrix @ samples_int_xim.T).T
            transformed_samples = [
                np.concatenate(
                    transform(
                        samples_int_xip[i],
                        samples_int_xim[i],
                        samples_rep_xip[i],
                        samples_rep_xim[i],
                    )
                )
                for i in tqdm.tqdm(range(n_samples_chunk), desc=f"Chunk {chunk_id}")
            ]
        else:
            transformed_samples = [
                np.concatenate(transform(samples_int_xip[i], samples_int_xim[i]))
                for i in tqdm.tqdm(range(n_samples_chunk), desc=f"Chunk {chunk_id}")
            ]
    else:
        samples_rep_xip = (binning_matrix @ samples_int_xip.T).T
        samples_rep_xim = (binning_matrix @ samples_int_xim.T).T

        transform = PureEBTransform(
            integration_method,
            theta=gg["meanr"],
            theta_int=gg_int["meanr"],
            tmin=min_sep,
            tmax=max_sep,
            # See mock_pure_eb_scatter: the single-grid operator evaluates at
            # nodes, so the pointwise operator chain snaps to the nearest fine
            # node in log theta. No effect on the dqags pointwise chain.
            snap_eval=integration_method == "fixed_operator",
            # The hybrid corner: the integral term is still the snapped
            # operator's, but the local 0.5*(xi_p +/- xi_m) term is read from
            # the reporting-grid draw rather than from one fine node.
            local_from_eval_grid=local_from_eval_grid,
        )
        if transform.snap_offset_max:
            print(
                f"[pure_eb] pointwise fixed_operator: evaluation points snapped "
                f"to the nearest fine node, max offset "
                f"{100 * transform.snap_offset_max:.3f}%",
                flush=True,
            )

        transformed_samples = [
            np.concatenate(
                transform(
                    samples_rep_xip[i],
                    samples_rep_xim[i],
                    samples_int_xip[i],
                    samples_int_xim[i],
                )
            )
            for i in tqdm.tqdm(range(n_samples_chunk), desc=f"Chunk {chunk_id}")
        ]

    eb_samples = np.array(transformed_samples)

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"pure_eb_chunk_{chunk_id}.npz")
    np.savez(out_path, eb_samples=eb_samples, chunk_id=chunk_id)
    print(f"Saved {n_samples_chunk} samples to {out_path}")
    return out_path


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chunk-id", type=int, required=True)
    ap.add_argument("--n-chunks", type=int, default=20)
    ap.add_argument("--n-samples", type=int, default=2000)
    ap.add_argument("--version", required=True)
    ap.add_argument("--blind", default="A")
    ap.add_argument("--cat-config", required=True)
    ap.add_argument("--xi-reporting", required=True)
    ap.add_argument("--xi-integration", required=True)
    ap.add_argument("--cov-integration", required=True)
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
    ap.add_argument(
        "--geometry",
        choices=("pointwise", "averaged"),
        default=None,
        help=(
            "Override the geometry the integration method defaults to. Omit "
            "for the mapping in pure_eb_transform.GEOMETRIES (dqags -> "
            "pointwise, fixed_operator -> averaged). 'averaged' with "
            "adaptive_dqags is the decision scan's cross term; 'pointwise' "
            "with fixed_operator is refused by the operator API."
        ),
    )
    ap.add_argument(
        "--local-from-eval-grid",
        action="store_true",
        help=(
            "Pointwise fixed_operator only: take the local 0.5*(xi_p +/- xi_m) "
            "term from the reporting-grid xi at the requested abscissa rather "
            "than from the snapped fine node, as the dqags pointwise chain "
            "does. The integral term still comes from the fine grid at the "
            "snapped node."
        ),
    )
    a = ap.parse_args(argv)
    compute_chunk(
        chunk_id=a.chunk_id,
        n_chunks=a.n_chunks,
        n_samples_total=a.n_samples,
        version=a.version,
        blind=a.blind,
        cat_config=a.cat_config,
        xi_reporting=a.xi_reporting,
        xi_integration=a.xi_integration,
        cov_integration=a.cov_integration,
        min_sep=a.min_sep,
        max_sep=a.max_sep,
        nbins=a.nbins,
        min_sep_int=a.min_sep_int,
        max_sep_int=a.max_sep_int,
        nbins_int=a.nbins_int,
        output_dir=a.out,
        geometry=a.geometry,
        integration_method=a.integration_method,
        rebinning_weights=a.rebinning_weights,
        integration_limits=a.integration_limits,
        local_from_eval_grid=a.local_from_eval_grid,
    )


if __name__ == "__main__":
    _from_cli()
