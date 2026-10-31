"""Recompute the calibrated min-PTE threshold on every scale-cut window.

The paper's min-PTE criterion is calibrated on the GLASS mocks at one angular
window: for each realisation take the minimum PTE across the six null-test
statistics, and read the 5% quantile of that distribution.  This script does
the same thing on every (theta_min, theta_max) pair of the 21-edge reporting
grid, with the COSEBIs *recomputed in each window* rather than held at their
fiducial values.

Per window:

* ``B_n<=6`` / ``B_n<=20`` -- COSEBIs recomputed from each mock's fine-binned
  xi_pm inside the window, with the COSEBI B covariance obtained by
  transforming that window's xi_pm covariance block through
  ``cosebis_covariance_from_xipm_covariance``.  No Hartlap factor, matching
  ``mock_campaign_analysis.py``.
* ``xi+^B`` / ``xi-^B`` / ``xi_tot^B`` -- the stored pure-mode B vectors on the
  20-bin reporting grid, cut by the bin-centre mask the pipeline uses, with
  the stored pure-mode covariance and the Hartlap factor that covariance
  artifact records (``pure_mode_covariance.hartlap_samples``: its Monte Carlo
  draw count, or none for the exact M C Mᵀ).
* ``C_ell^BB`` -- does not depend on the angular window; its per-mock PTE is
  taken unchanged from the per-mock statistics npz.

The COSEBI roots and norms are a function of theta_max/theta_min alone, and
across realisations that ratio varies only at the 1e-7 level, so each window
needs two root solves for the mock vectors (interpolated per realisation) and
one for the covariance rather than 351.  Everything else mirrors the pipeline
call for call; at the fiducial and full-range windows the result reproduces the
pipeline's own per-mock PTEs to machine precision.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from statistics import NormalDist

os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numba as nb  # noqa: E402
import numpy as np  # noqa: E402
import treecorr  # noqa: E402
from scipy import stats  # noqa: E402

from cosmo_numba.B_modes.cosebis import COSEBIS  # noqa: E402
from sp_validation.b_modes import scale_cut_to_bins  # noqa: E402

sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from pure_mode_covariance import hartlap_samples  # noqa: E402


STAT_NAMES = ("B_n<=6", "B_n<=20", "xi+^B", "xi-^B", "xi_tot^B", "C_ell^BB")
BOOTSTRAP_SAMPLES = 10_000
BOOTSTRAP_SEED = 42
CL_STAT_INDEX = 5

MIN_COMPARABLE_BINS = 2
"""Reporting bins a window must hold before its COSEBIs are comparable."""


# ── grid helpers ────────────────────────────────────────────────────────────


class _Grid:
    """The three attributes ``scale_cut_to_bins`` reads off a GGCorrelation."""

    def __init__(self, meanr, left_edges, right_edges):
        self.meanr = meanr
        self.left_edges = left_edges
        self.right_edges = right_edges


def _chi2_pte(data, covariance, n_samples=None):
    chi2 = float(data @ np.linalg.solve(covariance, data))
    dof = len(data)
    if n_samples is not None:
        chi2 *= (n_samples - dof - 2) / (n_samples - 1)
    return chi2, float(stats.chi2.sf(chi2, dof)), dof


def _threshold_with_interval(values, alpha):
    """The pipeline's bootstrap interval on the alpha-quantile."""

    threshold = float(np.quantile(values, alpha))
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, values.size, size=(BOOTSTRAP_SAMPLES, values.size))
    bootstrap = np.quantile(values[indices], alpha, axis=1)
    return threshold, [
        float(np.quantile(bootstrap, 0.16)),
        float(np.quantile(bootstrap, 0.84)),
    ]


def _order_statistic_interval(values, alpha):
    """Distribution-free 68% interval on the alpha-quantile.

    The number of sample points below the true alpha-quantile is
    Binomial(n, alpha); the 16th/84th percentiles of that count index the
    order statistics that bracket it.
    """

    ordered = np.sort(np.asarray(values, dtype=float))
    n = ordered.size
    lower = int(np.clip(stats.binom.ppf(0.16, n, alpha), 1, n)) - 1
    upper = int(np.clip(stats.binom.ppf(0.84, n, alpha) + 1, 1, n)) - 1
    return [float(ordered[lower]), float(ordered[upper])]


def _wilson_68(count, n):
    z = NormalDist().inv_cdf(0.84)
    fraction = count / n
    denominator = 1.0 + z * z / n
    centre = (fraction + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        fraction * (1.0 - fraction) / n + z * z / (4.0 * n * n)
    ) / denominator
    return [max(0.0, centre - half), min(1.0, centre + half)]


# ── the shared, read-once inputs ────────────────────────────────────────────


_STATE = {}


def load_inputs(args):
    """Everything the per-window work needs, read once."""

    with np.load(args.per_mock_statistics, allow_pickle=False) as data:
        stat_names = [str(name) for name in data["stat_names"]]
        pte_pipeline = np.asarray(data["pte"], dtype=float)
        pte_full_pipeline = np.asarray(data["pte_full"], dtype=float)
        cl_pte = np.asarray(data["pte"], dtype=float)[:, CL_STAT_INDEX]
        mock_ids = [str(value) for value in data["mock_ids"]]
        reporting_theta = np.asarray(data["theta"], dtype=float)
    if stat_names != list(STAT_NAMES):
        raise SystemExit(f"unexpected stat_names in {args.per_mock_statistics}")

    with np.load(args.mock_xi_fine, allow_pickle=False) as data:
        fine = {
            "mock_ids": [str(value) for value in data["mock_ids"]],
            "meanr": np.asarray(data["meanr"], dtype=float),
            "xip": np.asarray(data["xip"], dtype=float),
            "xim": np.asarray(data["xim"], dtype=float),
            "left_edges": np.asarray(data["left_edges"], dtype=float),
            "right_edges": np.asarray(data["right_edges"], dtype=float),
        }

    with np.load(args.mock_pure_eb_modes, allow_pickle=False) as data:
        pure = {
            "mock_ids": [str(value) for value in data["mock_ids"]],
            "theta": np.asarray(data["theta"], dtype=float),
            "xip_B": np.asarray(data["xip_B"], dtype=float),
            "xim_B": np.asarray(data["xim_B"], dtype=float),
        }

    # The mock ensembles carry their own ordering; the statistics npz is
    # ordered by sorted, zero-padded id.  Reindex both onto that order so a
    # row always means the same realisation.
    def reindex(block):
        lookup = {mock_id.zfill(5): i for i, mock_id in enumerate(block["mock_ids"])}
        missing = [m for m in mock_ids if m not in lookup]
        if missing:
            raise SystemExit(f"realisations missing from an ensemble: {missing[:5]}")
        return np.asarray([lookup[m] for m in mock_ids], dtype=int)

    fine_order = reindex(fine)
    pure_order = reindex(pure)

    with np.load(args.resolved_pure_covariance, allow_pickle=False) as data:
        if "cov_pure_eb" not in data:
            raise SystemExit(
                f"{args.resolved_pure_covariance} has no cov_pure_eb key"
            )
        pure_cov = np.asarray(data["cov_pure_eb"], dtype=float)
        pure_hartlap_n = hartlap_samples(data)
    nbins = pure["theta"].shape[1]
    if pure_cov.shape != (6 * nbins, 6 * nbins):
        raise SystemExit(f"unexpected pure-mode covariance shape {pure_cov.shape}")

    gg = treecorr.GGCorrelation(
        min_sep=args.xi_min_sep,
        max_sep=args.xi_max_sep,
        nbins=args.xi_nbins,
        sep_units="arcmin",
    )
    gg.read(str(args.cosebis_xi_grid))
    cov_xipm = np.asarray(np.loadtxt(args.cosebis_covariance), dtype=float)
    n_xi = len(gg.meanr)
    if cov_xipm.shape != (2 * n_xi, 2 * n_xi):
        raise SystemExit(
            f"expected a {2 * n_xi}x{2 * n_xi} xi covariance, got {cov_xipm.shape}"
        )
    if not np.allclose(gg.left_edges, fine["left_edges"], rtol=0, atol=1e-10):
        raise SystemExit(
            "the covariance xi grid and the mock fine grid do not share bin edges; "
            "the window would select different bins on either side"
        )

    return {
        "mock_ids": mock_ids,
        "pte_pipeline": pte_pipeline,
        "pte_full_pipeline": pte_full_pipeline,
        "cl_pte": cl_pte,
        "reporting_theta": reporting_theta,
        "fine_meanr": fine["meanr"][fine_order],
        "fine_xip": fine["xip"][fine_order],
        "fine_xim": fine["xim"][fine_order],
        "left_edges": fine["left_edges"],
        "right_edges": fine["right_edges"],
        "pure_theta": pure["theta"][pure_order],
        "pure_xip_B": pure["xip_B"][pure_order],
        "pure_xim_B": pure["xim_B"][pure_order],
        "pure_cov_xip": pure_cov[2 * nbins : 3 * nbins, 2 * nbins : 3 * nbins],
        "pure_cov_xim": pure_cov[3 * nbins : 4 * nbins, 3 * nbins : 4 * nbins],
        "pure_cov_cross": pure_cov[2 * nbins : 3 * nbins, 3 * nbins : 4 * nbins],
        "pure_hartlap_n": pure_hartlap_n,
        "xi_meanr": np.asarray(gg.meanr, dtype=float),
        "cov_xipm": cov_xipm,
        "nmodes": args.nmodes,
        "precision": args.precision,
        "cosebi_roots": args.cosebi_roots,
    }


def _roots_for_ratio(ratio, nmodes, precision):
    """COSEBI roots/norms for one theta_max/theta_min ratio.

    ``compute_roots_norms`` reads its two theta arguments only through their
    ratio, so anchoring at theta_min = 1 loses nothing.
    """

    obj = COSEBIS(1.0, float(ratio), nmodes, precision=precision)
    obj.compute_roots()
    return [np.asarray(entry, dtype=float) for entry in obj.roots], np.asarray(
        obj.norms, dtype=float
    )


# ── one window ──────────────────────────────────────────────────────────────


def cosebis_for_window(state, inds, theta_min_edge, theta_max_edge):
    """The two COSEBI statistics' per-mock PTEs inside one window.

    Mirrors ``mock_cosebis_scatter.compute_from_grid`` for the data vectors and
    ``mock_campaign_analysis._load_cosebis_covariance`` for the covariance --
    including the fact that the two use different theta grids (each mock's own
    pair-weighted meanr, and the data catalogue's).
    """

    nmodes = state["nmodes"]
    precision = state["precision"]
    theta_cov = state["xi_meanr"][inds]
    n_xi = state["xi_meanr"].size
    cov_inds = np.concatenate([inds, inds + n_xi])
    cov_obj = COSEBIS(
        float(theta_cov.min()), float(theta_cov.max()), nmodes, precision=precision
    )
    cov_cosebis = cov_obj.cosebis_covariance_from_xipm_covariance(
        theta_cov, state["cov_xipm"][cov_inds[:, None], cov_inds]
    )
    cov_B = cov_cosebis[nmodes:, nmodes:]

    theta_mocks = state["fine_meanr"][:, inds]
    xip = state["fine_xip"][:, inds]
    xim = state["fine_xim"][:, inds]
    n_mocks = theta_mocks.shape[0]
    tmins = theta_mocks.min(axis=1)
    tmaxs = theta_mocks.max(axis=1)

    # The pipeline gives every realisation its own COSEBIS object, built on
    # that realisation's own pair-weighted theta range.  Solving the roots 350
    # times per window is not affordable (~13 s each), but the roots and norms
    # are a function of theta_max/theta_min alone, and across realisations that
    # ratio spans a relative range of ~4e-7.  Solving at the two ends and
    # interpolating linearly in the ratio reproduces the pipeline's B_n to
    # ~1e-22 -- six orders of magnitude better than reusing one realisation's
    # roots, and far below the level at which any PTE moves.
    ratios = tmaxs / tmins
    low, high = float(ratios.min()), float(ratios.max())
    roots_low, norms_low = _roots_for_ratio(low, nmodes, precision)
    interpolate = state["cosebi_roots"] == "interpolate_ratio" and high > low
    if interpolate:
        roots_high, norms_high = _roots_for_ratio(high, nmodes, precision)

    mode_obj = COSEBIS(float(tmins[0]), float(tmaxs[0]), nmodes, precision=precision)
    bn = np.empty((n_mocks, nmodes))
    for index in range(n_mocks):
        if interpolate:
            weight = (float(ratios[index]) - low) / (high - low)
            mode_obj.roots = nb.typed.List(
                [a + (b - a) * weight for a, b in zip(roots_low, roots_high)]
            )
            mode_obj.norms = norms_low + (norms_high - norms_low) * weight
        else:
            mode_obj.roots = nb.typed.List(list(roots_low))
            mode_obj.norms = norms_low
        # The z = log(theta / tmin) mapping uses each realisation's own range,
        # exactly as its own COSEBIS object would have.
        mode_obj.tmin = float(tmins[index])
        mode_obj.tmax = float(tmaxs[index])
        bn[index] = mode_obj.cosebis_from_xipm(
            theta_mocks[index], xip[index], xim[index], parallel=True
        )[1]

    pte_6 = np.asarray(
        [_chi2_pte(row[:6], cov_B[:6, :6])[1] for row in bn], dtype=float
    )
    pte_20 = np.asarray([_chi2_pte(row, cov_B)[1] for row in bn], dtype=float)
    return pte_6, pte_20, {
        "n_xi_bins": int(inds.size),
        "cosebi_theta_min": float(theta_cov.min()),
        "cosebi_theta_max": float(theta_cov.max()),
    }


def pure_for_window(state, theta_min_edge, theta_max_edge):
    """The three pure-mode statistics' per-mock PTEs inside one window."""

    # The pipeline builds the bin-centre mask once, from the first
    # realisation's reporting theta, and applies it to every realisation.
    theta = state["pure_theta"][0]
    mask = (theta >= theta_min_edge) & (theta <= theta_max_edge)
    if not mask.any():
        return None, None, None, {"n_reporting_bins": 0}
    cov_xip = state["pure_cov_xip"][np.ix_(mask, mask)]
    cov_xim = state["pure_cov_xim"][np.ix_(mask, mask)]
    cov_cross = state["pure_cov_cross"][np.ix_(mask, mask)]
    cov_joint = np.block([[cov_xip, cov_cross], [cov_cross.T, cov_xim]])

    xip = state["pure_xip_B"][:, mask]
    xim = state["pure_xim_B"][:, mask]
    n_h = state["pure_hartlap_n"]
    pte_p = np.asarray(
        [_chi2_pte(row, cov_xip, n_h)[1] for row in xip], dtype=float
    )
    pte_m = np.asarray(
        [_chi2_pte(row, cov_xim, n_h)[1] for row in xim], dtype=float
    )
    pte_t = np.asarray(
        [
            _chi2_pte(np.concatenate([a, b]), cov_joint, n_h)[1]
            for a, b in zip(xip, xim)
        ],
        dtype=float,
    )
    return pte_p, pte_m, pte_t, {"n_reporting_bins": int(mask.sum())}


def window_ptes(state, edges, i, j):
    """The 350x6 PTE table for window (edges[i], edges[j])."""

    theta_min_edge, theta_max_edge = float(edges[i]), float(edges[j])
    grid = _Grid(state["fine_meanr"][0], state["left_edges"], state["right_edges"])
    start, stop = scale_cut_to_bins(grid, theta_min_edge, theta_max_edge)
    inds = np.arange(start, stop)
    info = {
        "edge_index_min": i,
        "edge_index_max": j,
        "theta_min": theta_min_edge,
        "theta_max": theta_max_edge,
    }
    n_mocks = state["fine_meanr"].shape[0]
    ptes = np.full((n_mocks, 6), np.nan)
    status = []

    if inds.size == 0:
        status.append("no fine xi bins in window")
    else:
        try:
            pte_6, pte_20, cosebi_info = cosebis_for_window(
                state, inds, theta_min_edge, theta_max_edge
            )
            ptes[:, 0], ptes[:, 1] = pte_6, pte_20
            info.update(cosebi_info)
        except Exception as error:  # noqa: BLE001 - recorded, not swallowed
            status.append(f"COSEBIs unavailable: {type(error).__name__}: {error}")
            info["n_xi_bins"] = int(inds.size)

    pte_p, pte_m, pte_t, pure_info = pure_for_window(
        state, theta_min_edge, theta_max_edge
    )
    info.update(pure_info)
    if 0 < info.get("n_reporting_bins", 0) < MIN_COMPARABLE_BINS:
        status.append(
            "one reporting bin: the COSEBI basis is solved over a single "
            "reporting bin's theta range, where it is too ill-conditioned to "
            "compare against the rest of the plane (median B_n<=20 PTE 0.13 "
            "here against ~0.45 on every window with two bins or more). "
            "Whether mpmath raises NoConvergence or returns these roots is "
            "environment-dependent, so the window is excluded structurally "
            "rather than by whether the solve happened to fail."
        )
    if pte_p is None:
        status.append("no reporting bins in window")
    else:
        ptes[:, 2], ptes[:, 3], ptes[:, 4] = pte_p, pte_m, pte_t

    ptes[:, 5] = state["cl_pte"]
    info["status"] = status
    info["statistics_available"] = [
        name for k, name in enumerate(STAT_NAMES) if np.all(np.isfinite(ptes[:, k]))
    ]
    return ptes, info


def summarise_window(ptes, info, alpha):
    available = np.asarray(
        [np.all(np.isfinite(ptes[:, k])) for k in range(len(STAT_NAMES))]
    )
    minima = np.min(ptes[:, available], axis=1)
    threshold, bootstrap = _threshold_with_interval(minima, alpha)
    return {
        **info,
        "n_statistics": int(available.sum()),
        "comparable": bool(
            available.all()
            and info.get("n_reporting_bins", 0) >= MIN_COMPARABLE_BINS
        ),
        "n_mocks": int(ptes.shape[0]),
        "threshold": threshold,
        "threshold_interval_68_bootstrap": bootstrap,
        "threshold_interval_68_order_statistic": _order_statistic_interval(
            minima, alpha
        ),
        "median_min_pte": float(np.median(minima)),
        # Exceedance fractions at the two fixed levels the response letter
        # quotes, alongside the calibrated quantile. Computed from the same
        # `minima` as `threshold` and `median_min_pte`: on a comparable window
        # this is the fraction of mocks whose minimum over all six statistics
        # falls below the level, and where a statistic is unavailable the
        # minimum is over the available ones, exactly as `threshold` already is.
        "alpha_0p01": float(np.mean(minima < 0.01)),
        "alpha_0p05": float(np.mean(minima < 0.05)),
        "median_pte": {
            name: (float(np.median(ptes[:, k])) if available[k] else None)
            for k, name in enumerate(STAT_NAMES)
        },
    }


# ── multiprocessing plumbing ────────────────────────────────────────────────


def _init_worker(state):
    _STATE["state"] = state


def _worker(task):
    state = _STATE["state"]
    edges, i, j, alpha = task
    ptes, info = window_ptes(state, edges, i, j)
    return summarise_window(ptes, info, alpha)


# ── validation ──────────────────────────────────────────────────────────────


def validate(state, edges, args):
    """Reproduce the pipeline's own per-mock PTEs at two known windows."""

    report = {}

    fid_ptes, fid_info = window_ptes(state, edges, *args.fiducial_window)
    reference = state["pte_pipeline"]
    report["fiducial"] = {
        "window_edges": list(args.fiducial_window),
        "theta_min": fid_info["theta_min"],
        "theta_max": fid_info["theta_max"],
        "reference": "mock_per_mock_statistics['pte']",
        "max_abs_difference": {
            name: float(np.max(np.abs(fid_ptes[:, k] - reference[:, k])))
            for k, name in enumerate(STAT_NAMES)
        },
        "status": fid_info["status"],
    }

    full_ptes, full_info = window_ptes(state, edges, 0, len(edges) - 1)
    reference_full = state["pte_full_pipeline"]
    report["full_range"] = {
        "window_edges": [0, len(edges) - 1],
        "theta_min": full_info["theta_min"],
        "theta_max": full_info["theta_max"],
        "reference": "mock_per_mock_statistics['pte_full']",
        "max_abs_difference": {
            name: float(np.max(np.abs(full_ptes[:, k] - reference_full[:, k])))
            for k, name in enumerate(STAT_NAMES)
        },
        "note": (
            "C_ell^BB is expected to disagree here by construction: this scan "
            "holds the harmonic statistic at its fixed 300<ell<1600 cut, while "
            "pte_full uses all 32 bands."
        ),
        "status": full_info["status"],
    }
    return report


# ── CLI ─────────────────────────────────────────────────────────────────────


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-mock-statistics", type=Path)
    parser.add_argument("--mock-xi-fine", type=Path)
    parser.add_argument("--mock-pure-eb-modes", type=Path)
    parser.add_argument(
        "--covariance",
        type=Path,
        help=(
            "xi_pm covariance on the fine grid, transformed to COSEBI space. "
            "The 'fiducial' reference covariance."
        ),
    )
    parser.add_argument("--pure-covariance", type=Path)
    parser.add_argument(
        "--covariance-mockmatched",
        type=Path,
        help="xi_pm covariance regenerated at the mocks' galaxy count.",
    )
    parser.add_argument(
        "--pure-covariance-mockmatched",
        type=Path,
        help="Pure-mode MC covariance regenerated at the mocks' galaxy count.",
    )
    parser.add_argument(
        "--reference-covariance",
        choices=("fiducial", "mockmatched"),
        default="fiducial",
        help=(
            "Which covariance pair the mock scatter is compared against, the "
            "same choice mock_campaign_analysis.py makes.  There is no C_l "
            "member here: the harmonic PTEs are read from the per-mock "
            "statistics npz, so they already carry that file's own choice."
        ),
    )
    parser.add_argument(
        "--cosebis-xi-grid",
        type=Path,
        help="TreeCorr fine-grid xi file supplying the covariance's theta axis.",
    )
    parser.add_argument("--xi-min-sep", type=float, default=0.5)
    parser.add_argument("--xi-max-sep", type=float, default=300.0)
    parser.add_argument("--xi-nbins", type=int, default=1000)
    parser.add_argument("--reporting-theta-min", type=float, default=1.0)
    parser.add_argument("--reporting-theta-max", type=float, default=250.0)
    parser.add_argument("--reporting-edges", type=int, default=21)
    parser.add_argument(
        "--fiducial-window",
        type=int,
        nargs=2,
        default=(9, 16),
        help="Edge indices of the fiducial window, for the validation gate.",
    )
    parser.add_argument("--global-alpha", type=float, default=0.05)
    parser.add_argument("--nmodes", type=int, default=20)
    parser.add_argument("--precision", type=int, default=120)
    parser.add_argument(
        "--cosebi-roots",
        choices=("interpolate_ratio", "shared"),
        default="interpolate_ratio",
        help=(
            "How each realisation's COSEBI roots/norms are obtained. "
            "'interpolate_ratio' solves at the two ends of the ensemble's "
            "theta_max/theta_min range and interpolates, reproducing the "
            "pipeline's per-realisation objects to ~1e-22 in B_n; 'shared' "
            "uses one solve for all realisations (~1e-16 in B_n, ~2e-5 in "
            "PTE). Each realisation keeps its own (theta_min, theta_max) "
            "either way."
        ),
    )
    parser.add_argument(
        "--emit",
        choices=("calibration", "figure", "validation"),
        default="calibration",
    )
    parser.add_argument(
        "--from-scan",
        type=Path,
        help=(
            "With --emit figure: read an already-written calibration json "
            "instead of recomputing the plane.  The figure is a pure function "
            "of that file, so the expensive pass runs once."
        ),
    )
    parser.add_argument(
        "--contour-level",
        type=float,
        default=0.01,
        help=(
            "Threshold level drawn as a labelled contour on the figure's first "
            "panel; pass a negative value to omit it."
        ),
    )
    parser.add_argument(
        "--highlight-theta-min",
        type=float,
        default=5.241,
        help="Lower-left corner of the region outlined on the figure.",
    )
    parser.add_argument(
        "--highlight-theta-max",
        type=float,
        default=36.196,
        help="Lower edge in theta_max of the region outlined on the figure.",
    )
    parser.add_argument("--out-file", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=1)
    return parser.parse_args(argv)


def _load_scan(path):
    """Read a calibration json back in, for the figure pass."""

    payload = json.loads(Path(path).read_text())
    if "windows" not in payload:
        raise SystemExit(f"{path} is not a window-scan calibration file")
    return payload


def _plot(payload, out_path, args):
    """The plane, with the two things that drive its extremes beside it.

    A lone heatmap of the 5% quantile invites the reader to take its extremes
    at face value.  They are not equally meaningful: the threshold is pushed
    *up* at short, narrow windows, where few and strongly correlated degrees of
    freedom keep every PTE near 0.5, and *down* at the widest windows, where
    the analytic COSEBI covariance understates the mock scatter and the B-mode
    PTEs are systematically low.  The two companion panels are exactly those
    two drivers, so the reader can see which part of the plane is calibration
    and which is artifact.
    """

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, TwoSlopeNorm

    try:
        from plotting_utils import PAPER_MPLSTYLE

        plt.style.use(PAPER_MPLSTYLE)
    except Exception as error:  # pragma: no cover - style is cosmetic
        print(f"Could not load the paper style ({error}); using defaults", flush=True)

    edges = np.asarray(payload["reporting_edges"], dtype=float)
    n = edges.size
    windows = payload["windows"]
    alpha = payload["global_alpha"]

    def plane(getter):
        grid = np.full((n - 1, n - 1), np.nan)
        for window in windows:
            i, j = window["edge_index_min"], window["edge_index_max"]
            value = getter(window)
            if value is not None:
                grid[i, j - 1] = value
        return grid

    usable = lambda w: w.get("comparable", w["n_statistics"] == len(STAT_NAMES))
    complete = [w for w in windows if usable(w)]
    threshold = plane(lambda w: w["threshold"] if usable(w) else None)
    dof = plane(lambda w: w["n_reporting_bins"])
    cosebi_pte = plane(lambda w: w["median_pte"]["B_n<=20"])
    missing = [
        (w["edge_index_min"], w["edge_index_max"])
        for w in windows
        if not usable(w)
    ]

    # pcolormesh is given explicit cell boundaries rather than centres, so the
    # hatched overlays land exactly on their cells; the contour is drawn on the
    # centres.
    x_centre, y_centre = edges[:-1], edges[1:]
    step = np.sqrt(edges[1] / edges[0])
    bounds = lambda c: np.concatenate(
        [[c[0] / step], np.sqrt(c[:-1] * c[1:]), [c[-1] * step]]
    )
    x_bound, y_bound = bounds(x_centre), bounds(y_centre)

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3))
    panels = [
        (
            threshold,
            dict(
                cmap="viridis",
                norm=LogNorm(
                    vmin=np.nanmin(threshold), vmax=np.nanmax(threshold)
                ),
            ),
            rf"calibrated threshold (${{\rm min}}$ PTE, $\alpha={alpha}$)",
        ),
        (
            dof,
            dict(cmap="cividis", vmin=1, vmax=np.nanmax(dof)),
            "reporting bins in window",
        ),
        (
            cosebi_pte,
            dict(cmap="RdBu_r", norm=TwoSlopeNorm(vmin=0.2, vcenter=0.5, vmax=0.6)),
            r"median mock PTE, $B_{n\leq 20}$",
        ),
    ]
    for ax, (grid, style, label) in zip(axes, panels):
        image = ax.pcolormesh(x_bound, y_bound, grid.T, shading="flat", **style)
        # A window holding one reporting bin gets a COSEBI basis solved over
        # that single bin's theta range, which is too ill-conditioned to
        # compare with the rest of the plane; hatch those rather than leave
        # holes, so they read as excluded rather than as missing data.
        for i_edge, j_edge in missing:
            ax.add_patch(
                plt.Rectangle(
                    (x_bound[i_edge], y_bound[j_edge - 1]),
                    x_bound[i_edge + 1] - x_bound[i_edge],
                    y_bound[j_edge] - y_bound[j_edge - 1],
                    facecolor="0.85",
                    edgecolor="0.15",
                    hatch="xxxx",
                    linewidth=0.4,
                    zorder=3,
                )
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\theta_{\rm min}$ [arcmin]")
        colorbar = fig.colorbar(image, ax=ax, pad=0.02)
        colorbar.set_label(label, fontsize=7)
        colorbar.ax.tick_params(labelsize=6)
        _mark_regions(ax, edges, payload, args)

    # The claim the calibration is quoted for is "threshold <= 0.01 over this
    # region"; drawing that level makes it checkable by eye instead of by
    # interpolating a logarithmic colourbar.
    level = args.contour_level
    if level is not None and np.nanmin(threshold) < level < np.nanmax(threshold):
        contour = axes[0].contour(
            x_centre,
            y_centre,
            threshold.T,
            levels=[level],
            colors="0.08",
            linewidths=1.1,
            zorder=6,
        )
        axes[0].clabel(
            contour, fmt={level: f"{level:g}"}, fontsize=6, inline=True
        )
    axes[0].set_ylabel(r"$\theta_{\rm max}$ [arcmin]")
    axes[0].legend(
        handles=[
            plt.Line2D(
                [], [], ls="none", marker="o", ms=4, mfc="none",
                mec="crimson", mew=1.2, label="fiducial window",
            ),
            plt.Line2D([], [], color="crimson", ls="--", lw=0.9, label="analysis region"),
            plt.Line2D([], [], color="0.08", lw=1.1, label=r"threshold $=0.01$"),
            plt.Rectangle(
                (0, 0), 1, 1, facecolor="0.85", edgecolor="0.15",
                hatch="xxxx", lw=0.4, label="single reporting bin",
            ),
        ],
        loc="lower right",
        fontsize=5.5,
        framealpha=0.9,
        handlelength=1.4,
    )
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(
        f"Saved {out_path} ({len(complete)} complete of {len(windows)} windows)",
        flush=True,
    )


def _mark_regions(ax, edges, payload, args):
    """Outline the region the analysis uses, and the fiducial window."""

    import matplotlib.pyplot as plt

    i, j = payload["fiducial_window_edges"]
    ax.plot(
        edges[i],
        edges[j],
        marker="o",
        ms=4,
        mfc="none",
        mec="crimson",
        mew=1.2,
        zorder=5,
    )
    lo, hi = args.highlight_theta_min, args.highlight_theta_max
    if lo is None or hi is None:
        return
    ax.add_patch(
        plt.Rectangle(
            (lo, hi),
            edges[-1] - lo,
            edges[-1] - hi,
            fill=False,
            edgecolor="crimson",
            linewidth=0.9,
            linestyle="--",
            zorder=4,
        )
    )


def main(argv=None):
    args = _parse_args(argv)
    if not 0.0 < args.global_alpha < 1.0:
        raise SystemExit("--global-alpha must lie strictly between 0 and 1")

    # The figure is a pure function of the calibration file, so this pass
    # opens none of the heavy inputs.
    if args.emit == "figure" and args.from_scan:
        _plot(_load_scan(args.from_scan), args.out_file, args)
        return

    # The declared `reference_covariance` decision picks the pair, exactly as
    # it does in mock_campaign_analysis.py.
    if args.reference_covariance == "mockmatched":
        args.cosebis_covariance = args.covariance_mockmatched
        args.resolved_pure_covariance = args.pure_covariance_mockmatched
    else:
        args.cosebis_covariance = args.covariance
        args.resolved_pure_covariance = args.pure_covariance
    required = {
        "--per-mock-statistics": args.per_mock_statistics,
        "--mock-xi-fine": args.mock_xi_fine,
        "--mock-pure-eb-modes": args.mock_pure_eb_modes,
        "--cosebis-xi-grid": args.cosebis_xi_grid,
        f"--covariance[-mockmatched] ({args.reference_covariance})":
            args.cosebis_covariance,
        f"--pure-covariance[-mockmatched] ({args.reference_covariance})":
            args.resolved_pure_covariance,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise SystemExit(
            "these inputs are required unless --emit figure --from-scan is "
            "used: " + ", ".join(missing)
        )
    state = load_inputs(args)
    edges = np.geomspace(
        args.reporting_theta_min, args.reporting_theta_max, args.reporting_edges
    )

    if args.emit == "validation":
        report = {
            "generated": datetime.now().isoformat(),
            "n_mocks": len(state["mock_ids"]),
            "cosebi_roots": args.cosebi_roots,
            **validate(state, edges, args),
        }
        args.out_file.parent.mkdir(parents=True, exist_ok=True)
        args.out_file.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return

    tasks = [
        (edges, i, j, args.global_alpha)
        for i in range(args.reporting_edges)
        for j in range(i + 1, args.reporting_edges)
    ]
    if args.processes > 1:
        import multiprocessing as mp

        with mp.get_context("fork").Pool(
            args.processes, initializer=_init_worker, initargs=(state,)
        ) as pool:
            windows = []
            for count, result in enumerate(pool.imap_unordered(_worker, tasks), 1):
                windows.append(result)
                print(
                    f"[{count}/{len(tasks)}] edges "
                    f"({result['edge_index_min']},{result['edge_index_max']}) "
                    f"threshold={result['threshold']:.5g}",
                    flush=True,
                )
    else:
        _init_worker(state)
        windows = []
        for count, task in enumerate(tasks, 1):
            windows.append(_worker(task))
            print(f"[{count}/{len(tasks)}] done", flush=True)
    windows.sort(key=lambda w: (w["edge_index_min"], w["edge_index_max"]))

    payload = {
        "generated": datetime.now().isoformat(),
        "stat_names": list(STAT_NAMES),
        "global_alpha": args.global_alpha,
        "nmodes": args.nmodes,
        "precision": args.precision,
        "cosebi_roots": args.cosebi_roots,
        "hartlap_n_pure": state["pure_hartlap_n"],
        "n_mocks": len(state["mock_ids"]),
        "reporting_edges": [float(value) for value in edges],
        "fiducial_window_edges": list(args.fiducial_window),
        "inputs": {
            "per_mock_statistics": str(args.per_mock_statistics),
            "mock_xi_fine": str(args.mock_xi_fine),
            "mock_pure_eb_modes": str(args.mock_pure_eb_modes),
            "reference_covariance": args.reference_covariance,
            "covariance": str(args.cosebis_covariance),
            "pure_covariance": str(args.resolved_pure_covariance),
            "cosebis_xi_grid": str(args.cosebis_xi_grid),
        },
        "windows": windows,
    }

    args.out_file.parent.mkdir(parents=True, exist_ok=True)
    if args.emit == "figure":
        _plot(payload, args.out_file, args)
        return
    args.out_file.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Saved {args.out_file}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
