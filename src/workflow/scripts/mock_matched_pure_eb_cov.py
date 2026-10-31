"""Build the mock-matched pure-E/B covariance.

Under ``--covariance-method exact_linear`` the covariance is ``M C_ξ Mᵀ`` of
the averaged transform (``pure_mode_covariance``), with C_ξ the mock-matched
ξ_± covariance; no draws are made. Under ``monte_carlo``, this is the paper's ``precompute_pure_eb_chunk.py`` plus
``gather_pure_eb_chunks.py`` path in one deterministic invocation.  The
processed 1000-bin ξ covariance is sampled in batches with seeds
``seed + batch_id`` — 20 batches of 100 under the defaults, which is
``precompute_pure_eb_chunk.py``'s own partition, so the two chains draw the
same streams; each sample is passed through Schneider (2022) under the
quadrature ``--integration-method`` selects, in the geometry that method
implies — rebinned to the 20-bin reporting grid then transformed pointwise
(``adaptive_dqags``), or transformed on the fine grid then averaged into the
reporting bins under the weights ``--rebinning-weights`` selects
(``fixed_operator``).  The output has the same keys and six-block
ordering as the paper covariance:

    [ξ+^E, ξ−^E, ξ+^B, ξ−^B, ξ+^amb, ξ−^amb]

The script adds only the sibling ``papers/bmodes/scripts`` directory to
``sys.path``, for the shared transform helper; it does not put the checkout's
``src/`` on the path.  The caller should provide the compatible cs_util
checkout via ``PYTHONPATH`` in the container invocation.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse

# The two script trees are separate sys.path roots; the transform helper that
# resolves `pure_eb_integration_method` lives with the paper-facing scripts and
# is shared rather than duplicated.
sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from pure_eb_transform import (  # noqa: E402
    AveragedPureEBTransform,
    PureEBTransform,
    add_integration_limits_argument,
    add_integration_method_argument,
    add_rebinning_weights_argument,
    check_rebinning_weights,
    geometry_for,
    resolve_integration_limits,
)
from pure_mode_covariance import (  # noqa: E402
    add_covariance_method_argument,
    covariance_record,
    exact_covariance,
)

_EB_KEYS = ("xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")


def _build_cosmology(cosmo_params):
    """Build the same CCL cosmology used by the paper chunk script."""

    import pyccl as ccl

    return ccl.Cosmology(
        Omega_c=cosmo_params["Omega_m"] - cosmo_params["Omega_b"],
        Omega_b=cosmo_params["Omega_b"],
        h=cosmo_params["h"],
        sigma8=cosmo_params["sigma_8"],
        n_s=cosmo_params["n_s"],
    )


def _load_xi(path, min_sep, max_sep, nbins, require_npairs=False):
    """Load a paper TreeCorr text dump and reconstruct its bin edges."""

    data = np.loadtxt(path, comments="#", max_rows=nbins)
    if data.shape[0] != nbins:
        raise ValueError(f"Expected {nbins} ξ rows in {path}, found {data.shape[0]}")
    if require_npairs and data.shape[1] <= 10:
        raise ValueError(f"Expected an npairs column at index 10 in {path}")
    edges = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins + 1)
    result = {
        "meanr": np.asarray(data[:, 1], dtype=float),
        "xip": np.asarray(data[:, 3], dtype=float),
        "xim": np.asarray(data[:, 4], dtype=float),
        "left_edges": edges[:-1],
        "right_edges": edges[1:],
    }
    if data.shape[1] > 10:
        result["npairs"] = np.asarray(data[:, 10], dtype=float)
    return result


def _load_nz(version, blind, cat_config):
    """Resolve the (version, blind) n(z) the way the paper chunk script does.

    ``precompute_pure_eb_chunk.py`` builds a ``CosmologyValidation`` over the
    single version and calls ``cv.get_redshift(version)``, which strips
    ``_leak_corr``, opens cat_config's ``redshift_path`` for the base version
    and substitutes the blind suffix.  Going through the same call is what
    keeps this chain's mean vector identical to the semi-analytic one: the two
    differ only in the covariance they sample.  ``CosmologyValidation`` insists
    on an output directory it can create, and this script has no product to put
    there, so it gets a throwaway one.
    """
    from sp_validation.cosmo_val import CosmologyValidation

    cv = CosmologyValidation(
        versions=[version],
        catalog_config=cat_config,
        output_dir=tempfile.mkdtemp(prefix="mock_matched_pure_eb_cov_"),
    )
    cv.blind = blind
    z, nz = cv.get_redshift(version)
    return np.column_stack([z, nz])


def _pointwise_binning_matrix(left_edges, right_edges, theta_int, nbins_int):
    """The published pointwise rebin of the *input* ξ_±: an unweighted mean.

    Verbatim the construction in ``precompute_pure_eb_chunk.py`` (lines
    129-146): digitize the fine grid into the reporting bins and average the
    members of each bin with equal weight.  This is a different operation from
    the fine-to-reporting average of the *pure modes* that
    ``pure_eb_rebinning_weights`` governs: it produces the reporting-grid ξ_±
    the pointwise transform is evaluated against, before any decomposition, and
    is what the published chain did.  It is deliberately not routed through
    ``pure_eb_transform.make_binning_matrix`` — under the pointwise geometry
    that helper raises, since ``none`` is the only weighting the geometry
    admits, and reaching for its averaged form here would make this covariance
    disagree with ``pure_eb_semianalytic_data`` for reasons that have nothing to
    do with the input covariance.
    """

    n_reporting = len(left_edges)
    reporting_bin_edges = np.concatenate([left_edges, [right_edges[-1]]])
    bin_indices = np.digitize(theta_int, reporting_bin_edges) - 1
    valid_mask = (bin_indices >= 0) & (bin_indices < n_reporting)
    row_indices, col_indices = bin_indices[valid_mask], np.where(valid_mask)[0]

    matrix = sparse.csr_matrix(
        (np.ones(len(row_indices)), (row_indices, col_indices)),
        shape=(n_reporting, nbins_int),
    )
    row_sums = np.array(matrix.sum(axis=1)).flatten()
    return sparse.diags(1 / row_sums) @ matrix


def _print_drop_counts(label, counts):
    details = ", ".join(f"{name}={count}" for name, count in zip(_EB_KEYS, counts))
    print(f"Dropped non-finite {label}: {details}", flush=True)


def _sample_bounds(n_samples, n_batches, batch_id):
    samples_per_batch = n_samples // n_batches
    start = batch_id * samples_per_batch
    stop = start + samples_per_batch if batch_id < n_batches - 1 else n_samples
    return start, stop


def build_covariance(
    *,
    covariance_path,
    xi_reporting_path,
    xi_integration_path,
    version,
    blind,
    cat_config,
    output_path,
    integration_method,
    rebinning_weights,
    covariance_method,
    n_samples=2000,
    n_batches=20,
    seed=42,
    integration_limits: str,
    min_sep=1.0,
    max_sep=250.0,
    nbins=20,
    min_sep_int=0.5,
    max_sep_int=300.0,
    nbins_int=1000,
    pad_xim=True,
    interp_order=5,
    diagonal_covariance=False,
    transform_averaged=False,
    geometry=None,
    local_from_eval_grid=False,
):
    """Generate and save the 120×120 mock-matched covariance."""

    if n_samples < 2 or n_batches < 1 or n_batches > n_samples:
        raise ValueError("Require 2 <= n_samples and 1 <= n_batches <= n_samples")

    # Resolved and checked up front, before the covariance load and the n(z)
    # lookup: a universe that pairs the quadrature and the weighting
    # incoherently should say so immediately.
    # See mock_pure_eb_scatter: resolved through geometry_for, not an OR, so
    # an explicit "pointwise" survives a fixed_operator method.
    geometry = geometry_for(
        integration_method, "averaged" if transform_averaged else geometry
    )
    averaged = geometry == "averaged"
    theta_window = resolve_integration_limits(integration_limits, geometry)
    check_rebinning_weights(rebinning_weights, geometry)

    from sp_validation.cosmology import PLANCK18, get_theo_xi

    reporting = _load_xi(xi_reporting_path, min_sep, max_sep, nbins)
    integration = _load_xi(
        xi_integration_path,
        min_sep_int,
        max_sep_int,
        nbins_int,
        require_npairs=True,
    )
    covariance = np.asarray(np.loadtxt(covariance_path), dtype=float)
    expected_shape = (2 * nbins_int, 2 * nbins_int)
    if covariance.shape != expected_shape:
        raise ValueError(
            f"Expected processed ξ covariance shape {expected_shape}, "
            f"found {covariance.shape}"
        )
    if diagonal_covariance:
        # Diagnostic only: retain the CosmoCov diagonal while removing all
        # fine-bin and xi+/xi- correlations.  The default is the full stored
        # covariance used for the paper product.
        covariance = np.diag(np.diag(covariance))

    z_dist = _load_nz(version, blind, cat_config)
    cosmo = _build_cosmology(dict(PLANCK18))
    mean_int = np.concatenate(
        get_theo_xi(
            theta=integration["meanr"],
            z=z_dist[:, 0],
            nz=z_dist[:, 1],
            backend="ccl",
            cosmo=cosmo,
        )
    )
    all_samples = []

    # One geometry for the whole run, so under `fixed_operator` the operator is
    # built once and reused across the 2000 draws and the data vector.
    # The pointwise chain always needs the uniform fine-to-reporting matrix;
    # so does the averaged reporting-local-term variant, which feeds the
    # transform its reporting draws.
    if not averaged or local_from_eval_grid:
        binning_matrix = _pointwise_binning_matrix(
            reporting["left_edges"],
            reporting["right_edges"],
            integration["meanr"],
            nbins_int,
        )
    if averaged:
        transform = AveragedPureEBTransform(
            integration_method,
            left_edges=reporting["left_edges"],
            right_edges=reporting["right_edges"],
            theta_int=integration["meanr"],
            rebinning_weights=rebinning_weights,
            npairs=integration.get("npairs"),
            theta_window=theta_window,
            pad_xim=pad_xim,
            interp_order=interp_order,
            parallel=True,
            local_from_eval_grid=local_from_eval_grid,
        )
        dropped_mc = np.zeros(len(_EB_KEYS), dtype=int)
    else:
        transform = PureEBTransform(
            integration_method,
            theta=reporting["meanr"],
            theta_int=integration["meanr"],
            tmin=min_sep,
            tmax=max_sep,
            pad_xim=pad_xim,
            interp_order=interp_order,
            parallel=True,
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
                "[pure_eb] pointwise fixed_operator: evaluation points snapped "
                "to the nearest fine node, max offset "
                f"{100 * transform.snap_offset_max:.3f}%",
                flush=True,
            )

    if covariance_method == "exact_linear":
        if not averaged or local_from_eval_grid:
            raise ValueError(
                "exact_linear needs the averaged transform on the fine grid alone "
                "(fixed_operator, no --local-from-eval-grid)"
            )
    # The exact covariance makes no draws; the transform above serves both.
    mc_batches = range(n_batches) if covariance_method == "monte_carlo" else ()
    for batch_id in mc_batches:
        start, stop = _sample_bounds(n_samples, n_batches, batch_id)
        n_batch = stop - start
        print(
            f"Batch {batch_id + 1}/{n_batches}: samples {start}:{stop} "
            f"({n_batch}), seed={seed + batch_id}",
            flush=True,
        )
        rng = np.random.default_rng(seed=seed + batch_id)
        samples_int = rng.multivariate_normal(mean_int, covariance, size=n_batch)
        samples_xip_int = samples_int[:, :nbins_int]
        samples_xim_int = samples_int[:, nbins_int:]
        if not averaged or local_from_eval_grid:
            # See precompute_pure_eb_chunk: under `uniform` weights these are
            # the same numbers the averaged transform's own local term would
            # be, so the substitution is a no-op on the MC draws; it bites on
            # the measured data vector below.
            samples_xip_rep = (binning_matrix @ samples_xip_int.T).T
            samples_xim_rep = (binning_matrix @ samples_xim_int.T).T

        transformed = []
        for index in range(n_batch):
            if averaged:
                modes = transform(
                    samples_xip_int[index],
                    samples_xim_int[index],
                    *(
                        (samples_xip_rep[index], samples_xim_rep[index])
                        if local_from_eval_grid
                        else ()
                    ),
                )
                dropped_mc += transform.last_dropped
            else:
                modes = transform(
                    samples_xip_rep[index],
                    samples_xim_rep[index],
                    samples_xip_int[index],
                    samples_xim_int[index],
                )
            transformed.append(np.concatenate(modes))
        all_samples.append(np.asarray(transformed, dtype=float))

    if covariance_method == "exact_linear":
        cov_pure_eb = exact_covariance(transform, covariance)
        record = covariance_record("exact_linear")
    else:
        if averaged:
            _print_drop_counts("MC outputs", dropped_mc)
        eb_samples = np.vstack(all_samples)
        cov_pure_eb = np.cov(eb_samples.T)
        record = covariance_record("monte_carlo", n_samples=n_samples)
    if cov_pure_eb.shape != (120, 120):
        raise ValueError(
            f"Expected 120x120 output covariance, found {cov_pure_eb.shape}"
        )

    # Match gather_pure_eb_chunks.py's package exactly, with provenance extras.
    if averaged:
        modes_data = transform(
            integration["xip"],
            integration["xim"],
            *(
                (reporting["xip"], reporting["xim"])
                if local_from_eval_grid
                else ()
            ),
        )
        _print_drop_counts("data outputs", transform.last_dropped)
    else:
        modes_data = transform(
            reporting["xip"],
            reporting["xim"],
            integration["xip"],
            integration["xim"],
        )
    xip_E, xim_E, xip_B, xim_B, xip_amb, xim_amb = modes_data
    if averaged and not all(
        np.all(np.isfinite(values)) for values in modes_data
    ):
        raise ValueError("Non-finite pure E/B output after averaging")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        theta=reporting["meanr"],
        theta_int=integration["meanr"],
        xip_total=reporting["xip"],
        xim_total=reporting["xim"],
        xip_E=xip_E,
        xim_E=xim_E,
        xip_B=xip_B,
        xim_B=xim_B,
        xip_amb=xip_amb,
        xim_amb=xim_amb,
        cov_pure_eb=cov_pure_eb,
        **record,
        **({"n_batches": n_batches, "seed": seed} if covariance_method == "monte_carlo" else {}),
        covariance_path=str(covariance_path),
        version=version,
        blind=blind,
        pad_xim=pad_xim,
        interp_order=interp_order,
        diagonal_covariance=diagonal_covariance,
        integration_limits=integration_limits,
        transform="averaged" if averaged else "pointwise",
        integration_method=integration_method,
        weighting=rebinning_weights,
    )
    print(f"Saved {output_path} with cov_pure_eb shape {cov_pure_eb.shape}", flush=True)
    return output_path


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cov-integration", required=True)
    parser.add_argument("--xi-reporting", required=True)
    parser.add_argument("--xi-integration", required=True)
    parser.add_argument(
        "--version",
        required=True,
        help="Catalog version whose n(z) cat_config resolves, e.g. "
        "SP_v1.4.6.3_leak_corr.",
    )
    parser.add_argument("--blind", default="A")
    parser.add_argument(
        "--cat-config", required=True, help="Path to cat_config.yaml"
    )
    parser.add_argument("--out", required=True, help="Final .npz output path")
    parser.add_argument("--n-samples", type=int, default=2000)
    parser.add_argument(
        "--n-batches",
        "--n-chunks",
        dest="n_batches",
        type=int,
        default=20,
        help=(
            "Deterministic RNG batches, not separate output files. The default "
            "matches precompute_pure_eb_chunk's 20 chunks: the batch id is the "
            "seed offset, so any other value draws a different RNG stream than "
            "the semi-analytic chain this one is meant to differ from only in "
            "its input covariance."
        ),
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-sep", type=float, default=1.0)
    parser.add_argument("--max-sep", type=float, default=250.0)
    parser.add_argument("--nbins", type=int, default=20)
    parser.add_argument("--min-sep-int", type=float, default=0.5)
    parser.add_argument("--max-sep-int", type=float, default=300.0)
    parser.add_argument("--nbins-int", type=int, default=1000)
    parser.add_argument(
        "--no-pad-xim",
        action="store_true",
        help="Diagnostic: disable cosmo_numba's default xi- padding.",
    )
    parser.add_argument(
        "--interp-order",
        type=int,
        default=5,
        choices=(1, 3, 5, 7, 9),
        help="Interpolation order passed to get_pure_EB_modes.",
    )
    parser.add_argument(
        "--diagonal-covariance",
        action="store_true",
        help="Diagnostic: zero all off-diagonal CosmoCov elements.",
    )
    parser.add_argument(
        "--geometry",
        choices=("pointwise", "averaged"),
        default=None,
        help="Override the geometry the integration method defaults to.",
    )
    parser.add_argument(
        "--transform-averaged",
        action="store_true",
        help=(
            "Force the averaged geometry under adaptive_dqags; fixed_operator "
            "already implies it."
        ),
    )
    parser.add_argument(
        "--local-from-eval-grid",
        action="store_true",
        help='Pointwise fixed_operator only: take the local 0.5*(xi_p +/- xi_m) term from the reporting-grid xi at the requested abscissa rather than from the snapped fine node, as the dqags pointwise chain does. The integral term still comes from the fine grid at the snapped node.',
    )
    add_integration_method_argument(parser)
    add_integration_limits_argument(parser)
    add_rebinning_weights_argument(parser)
    add_covariance_method_argument(parser)
    return parser.parse_args(argv)


def _main(argv=None):
    args = _parser(argv)
    build_covariance(
        covariance_path=args.cov_integration,
        xi_reporting_path=args.xi_reporting,
        xi_integration_path=args.xi_integration,
        version=args.version,
        blind=args.blind,
        cat_config=args.cat_config,
        output_path=args.out,
        integration_method=args.integration_method,
        rebinning_weights=args.rebinning_weights,
        covariance_method=args.covariance_method,
        n_samples=args.n_samples,
        n_batches=args.n_batches,
        seed=args.seed,
        integration_limits=args.integration_limits,
        min_sep=args.min_sep,
        max_sep=args.max_sep,
        nbins=args.nbins,
        min_sep_int=args.min_sep_int,
        max_sep_int=args.max_sep_int,
        nbins_int=args.nbins_int,
        pad_xim=not args.no_pad_xim,
        interp_order=args.interp_order,
        diagonal_covariance=args.diagonal_covariance,
        transform_averaged=args.transform_averaged,
        geometry=args.geometry,
        local_from_eval_grid=args.local_from_eval_grid,
    )


if __name__ == "__main__":
    _main()
