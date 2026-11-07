"""Analyse the six Paper-II B-mode statistics across GLASS mocks.

It reads the mock ensembles produced by the scatter sweeps, compares them with
the paper's stored analytic covariances, and reports the result as one metric
file and three figures.

The three covariance conventions are intentionally kept separate:

* COSEBIs: the paper 0.5--300 arcmin, 1000-bin covariance is transformed with
  the same ``COSEBIS.cosebis_covariance_from_xipm_covariance`` call used by the
  paper scripts; no Hartlap factor is applied.
* Pure ξ E/B: the stored 120 x 120 six-block covariance is sliced to the
  [ξ+^B, ξ-^B] blocks and the 12--83 arcmin cuts; the Hartlap factor is the
  one the covariance artifact records (its Monte Carlo draw count, or none for
  the exact M C Mᵀ; ``pure_mode_covariance.hartlap_samples``).
* C_l^BB: the stored 32 x 32 covariance is selected with the same powspace
  bin-edge logic as the paper scripts; no Hartlap factor is applied.  Where a
  BB covariance computed directly on a mock catalogue exists it is passed on
  ``--cl-cov-mockmatched`` and used as-is; otherwise ``--cl-cov-scale``
  rescales the stored data covariance by the squared data/mock galaxy-count
  ratio, which is exact for BB because these mocks carry no B signal.

Compute and reporting are separate invocations, because they are separate
artifacts:

* ``--emit-stats PATH`` is the compute pass.  It reads the mock ensembles and
  the reference covariances, forms the six statistics for every realisation,
  and writes *one* npz — the per-realisation table plus the analytic σ vectors
  and cut masks the reporting passes need.  This is the expensive step and the
  only one that touches a covariance.
* ``--from-stats PATH --emit {evidence,sigma-ratio-figure,pte-uniformity-figure,
  correlation-figure}`` are the cheap derivations.  Each reads that npz and
  writes exactly one file.  They recompute nothing: every number they report is
  a function of the stored arrays.

``--reference-covariance`` selects the covariance triple through the declared
``reference_covariance`` decision: ``data_covariance`` uses the covariances passed on
the base flags, ``mockmatched`` the ``*-mockmatched`` companions.  Per-mock
products arrive as stacked ensembles.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import treecorr
from astropy.io import fits
from cosmo_numba.B_modes.cosebis import COSEBIS
from scipy import stats

from sp_validation.b_modes import scale_cut_to_bins

sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from pure_mode_covariance import MC_SAMPLES_WITHOUT_RECORD, hartlap_samples  # noqa: E402

from numerics import (  # noqa: E402
    compute_chi2_pte,
    ell_bin_mask,
    get_powspace_bin_edges,
)
from plotting_utils import (  # noqa: E402
    CATEGORICAL_COLORS,
    PAPER_MPLSTYLE,
    SCATTER_ERRORBAR,
)

plt.style.use(PAPER_MPLSTYLE)


STAT_NAMES = (
    "B_n<=6",
    "B_n<=20",
    "xi+^B",
    "xi-^B",
    "xi_tot^B",
    "C_ell^BB",
)
COSEBI_THETA_MIN = 250.0 ** (9.0 / 20.0)
COSEBI_THETA_MAX = 250.0 ** (16.0 / 20.0)
FULL_COSEBI_THETA_MIN = 1.0
FULL_COSEBI_THETA_MAX = 250.0
# The two readings of the paper's fiducial θ cut, as the root
# `scale_cut_convention` decision spells them: the grid bin edges the pipeline
# actually cut at, and the rounded arcmin values the paper text prints.
COSEBI_THETA_BY_CONVENTION = {
    "grid_bin_edges": (COSEBI_THETA_MIN, COSEBI_THETA_MAX),
    "rounded_arcmin": (12.0, 83.0),
}
# The pure-mode blocks live on the 20-bin reporting grid and are cut by a
# bin-CENTRE mask, so both conventions select the same seven bins (centres
# 13.99 .. 73.29 arcmin; the nearest excluded centres are 10.61 and 96.58, far
# from either boundary pair).  The cut is still resolved from the declared
# convention rather than pinned here, so that the spec stays the single source
# of truth and a future change to the reporting grid cannot silently desync it.
PURE_THETA_BY_CONVENTION = dict(COSEBI_THETA_BY_CONVENTION)
PURE_THETA_MIN, PURE_THETA_MAX = PURE_THETA_BY_CONVENTION["grid_bin_edges"]
ELL_MIN = 300.0
ELL_MAX = 1600.0
# These two are galaxy COUNTS in disguise, not measured sky densities.  Both are
# a catalogue's galaxy count divided by the same fiducial 2894 deg^2 area
# convention: 4.9573 arcmin^-2 is the data's 51.65M weight-effective galaxies and
# 3.3697 arcmin^-2 is the v1.4.6 GLASS suite's 35 078 729 unit-weight galaxies.
# The mock and data catalogues occupy the same sky, so the shared area cancels
# and only their ratio is physical: it is the count ratio 51.65M / 35.08M =
# 1.4711, i.e. the v1.4.6 suite is that much NOISIER than the data.  Neither
# number is the density at which GLASS was asked to sample (the v1.4.6 suite
# covers ~2.6e3 deg^2 of sky, so its density there is ~3.7 arcmin^-2), and
# neither applies to the v1.4.6.3_v2 suite, whose 60 423 575 galaxies make it
# QUIETER than the data by the same kind of ratio (51.65/60.42 = 0.855).  The v2
# suite's C_l^BB reference is a directly computed mock covariance passed on
# --cl-cov-mockmatched, not a rescale of the data's.
FIDUCIAL_NEFF = 4.957279270321334
MOCK_NEFF = 3.369700396483938
MOCKMATCHED_CL_COV_SCALE = (FIDUCIAL_NEFF / MOCK_NEFF) ** 2
STAT_TEX = {
    "B_n<=6": r"$B_{n\leq 6}$",
    "B_n<=20": r"$B_{n\leq 20}$",
    "xi+^B": r"$\xi_+^{B}$",
    "xi-^B": r"$\xi_-^{B}$",
    "xi_tot^B": r"$\xi_{\rm tot}^{B}$",
    "C_ell^BB": r"$C_\ell^{BB}$",
}


def _native(values):
    return np.asarray(values, dtype=float)


def _normalise_ids(ids):
    return sorted({str(mock_id).zfill(5) for mock_id in ids})


def _default_paths():
    repo = "/automnt/n17data/cdaley/unions/code/sp_validation"
    analysis = "results/glass_mock"
    cov_base = f"{repo}/cosmo_inference/data/covariance"
    cosebis_cov_base = (
        "covariance_SP_v1.4.6.3_leak_corr_A_g_minsep=0.5_maxsep=300.0_nbins=1000_masked"
    )
    return {
        "cl_dir": "/n09data/guerrini/glass_mock_v1.4.6/results",
        "pseudo_cl_dir": analysis,
        "cosebis_xi_grid": (
            "/n17data/cdaley/unions/code/sp_validation/cosmo_val/output/"
            "SP_v1.4.6.3_leak_corr_xi_minsep=0.5_maxsep=300.0_"
            "nbins=1000_npatch=1.txt"
        ),
        "cosebis_cov": f"{cov_base}/{cosebis_cov_base}/{cosebis_cov_base}_processed.txt",
        "pure_cov": (
            "/automnt/n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/"
            "results/paper_plots/intermediate/"
            "SP_v1.4.6.3_leak_corr_A_pure_eb_semianalytic.npz"
        ),
        "cl_cov": (
            "/n17data/cdaley/unions/code/sp_validation/cosmo_val/output/"
            "pseudo_cl_cov_SP_v1.4.6.3_leak_corr_blind=A_powspace_nbins=32.fits"
        ),
    }


def _path_for(kind: str, directory: str | Path, mock_id: str) -> Path:
    directory = Path(directory)
    patterns = {
        "cosebis": f"cosebis_glass_mock_{mock_id}.npz",
        "pure": f"pure_eb_glass_mock_{mock_id}.npz",
        "fine": f"gg_glass_mock_{mock_id}_nbins=1000.fits",
        "cl": f"cl_glass_mock_{mock_id}_4096.npy",
        "pseudo_cl": f"pseudo_cl_glass_mock_{mock_id}_powspace_nbins=32.fits",
    }
    return directory / patterns[kind]


def _read_ensemble(path):
    """Load a stacked mock ensemble npz into a plain dict of arrays."""

    with np.load(path, allow_pickle=False) as data:
        ensemble = {key: data[key] for key in data.files}
    if "mock_ids" not in ensemble:
        raise KeyError(f"{path} is not a mock ensemble; keys={sorted(ensemble)}")
    ensemble["mock_ids"] = [str(value) for value in np.asarray(ensemble["mock_ids"])]
    ensemble["source"] = str(path)
    return ensemble


def _load_cosebis(ensemble, index):
    source = f"{ensemble['source']}[{ensemble['mock_ids'][index]}]"
    if "Bn" not in ensemble:
        raise KeyError(f"{source} has no Bn key")
    bn = _native(ensemble["Bn"][index])
    en = _native(ensemble["En"][index]) if "En" in ensemble else None
    theta_min = (
        float(np.asarray(ensemble["theta_min_actual"][index]).squeeze())
        if "theta_min_actual" in ensemble
        else None
    )
    theta_max = (
        float(np.asarray(ensemble["theta_max_actual"][index]).squeeze())
        if "theta_max_actual" in ensemble
        else None
    )
    if bn.shape != (20,):
        raise ValueError(f"Expected 20 COSEBI B modes in {source}, found {bn.shape}")
    return {
        "Bn": bn,
        "En": en,
        "stored_theta_min": theta_min,
        "stored_theta_max": theta_max,
        "source": source,
    }


def _stored_cosebis_is_fiducial(cosebis):
    """Allow bin-centre rounding, but reject the current 1--250 products."""

    tmin, tmax = cosebis["stored_theta_min"], cosebis["stored_theta_max"]
    if tmin is None or tmax is None:
        return False
    return abs(tmin - COSEBI_THETA_MIN) < 1.0 and abs(tmax - COSEBI_THETA_MAX) < 2.0


class _Grid:
    """The three attributes `scale_cut_to_bins` reads off a GGCorrelation."""

    def __init__(self, meanr, left_edges, right_edges):
        self.meanr = meanr
        self.left_edges = left_edges
        self.right_edges = right_edges


def _recompute_cosebis_from_fine(fine_ensemble, index):
    """Recompute stale full-range mock COSEBIs on the paper fiducial cut."""

    source = f"{fine_ensemble['source']}[{fine_ensemble['mock_ids'][index]}]"
    grid = _Grid(
        fine_ensemble["meanr"][index],
        fine_ensemble["left_edges"],
        fine_ensemble["right_edges"],
    )
    start, stop = scale_cut_to_bins(grid, COSEBI_THETA_MIN, COSEBI_THETA_MAX)
    inds = np.arange(start, stop)
    theta = _native(grid.meanr[inds])
    xip = _native(fine_ensemble["xip"][index][inds])
    xim = _native(fine_ensemble["xim"][index][inds])
    if len(theta) == 0:
        raise ValueError(f"No fine ξ bins in the COSEBI fiducial cut for {source}")
    cosebis = COSEBIS(float(theta.min()), float(theta.max()), 20, precision=120)
    en, bn = cosebis.cosebis_from_xipm(theta, xip, xim, parallel=True)
    return {
        "Bn": _native(bn),
        "En": _native(en),
        "stored_theta_min": float(theta.min()),
        "stored_theta_max": float(theta.max()),
        "source": source,
    }


def _load_or_recompute_cosebis(cosebis_ensemble, fine_ensemble, index):
    stored = _load_cosebis(cosebis_ensemble, index)
    if _stored_cosebis_is_fiducial(stored):
        stored["scale_action"] = "used_stored_fiducial_vector"
        return stored
    if fine_ensemble is None:
        raise ValueError(
            f"{stored['source']} stores a non-fiducial/unknown COSEBI cut "
            f"({stored['stored_theta_min']}, {stored['stored_theta_max']}); "
            "no fine ξ ensemble was given for recomputation"
        )
    print(
        f"COSEBI cut mismatch for {stored['source']}: stored "
        f"[{stored['stored_theta_min']}, {stored['stored_theta_max']}], "
        f"recomputing at [{COSEBI_THETA_MIN:.6g}, {COSEBI_THETA_MAX:.6g}]",
        flush=True,
    )
    result = _recompute_cosebis_from_fine(fine_ensemble, index)
    result["stored_theta_min"] = stored["stored_theta_min"]
    result["stored_theta_max"] = stored["stored_theta_max"]
    result["scale_action"] = "recomputed_from_fine_xi"
    return result


def _load_pure(ensemble, index):
    source = f"{ensemble['source']}[{ensemble['mock_ids'][index]}]"
    required = ("theta", "xip_B", "xim_B")
    missing = [key for key in required if key not in ensemble]
    if missing:
        raise KeyError(f"{source} is missing pure-EB keys {missing}")
    result = {key: _native(ensemble[key][index]) for key in required}
    # The reporting binning and the integration cut are ensemble-wide, not
    # per-realisation: they are stored once by mock_pure_eb_scatter.
    for key in ("reporting_left_edges", "reporting_right_edges"):
        if key in ensemble:
            result[key] = _native(ensemble[key])
    if "integration_limits" in ensemble:
        result["integration_limits"] = str(
            np.asarray(ensemble["integration_limits"]).squeeze()
        )
    if result["theta"].shape != (20,):
        raise ValueError(f"Expected 20 pure-EB reporting bins in {source}")
    numeric_keys = (
        "theta",
        "xip_B",
        "xim_B",
        "reporting_left_edges",
        "reporting_right_edges",
    )
    if not all(np.all(np.isfinite(result[key])) for key in numeric_keys if key in result):
        raise ValueError(f"Non-finite pure-EB vector in {source}")
    return result


def _load_cl(path: Path):
    array = _native(np.load(path, allow_pickle=False))
    if array.ndim != 2 or array.shape[0] < 5 or array.shape[1] != 32:
        raise ValueError(
            f"Expected external C_l array with at least 5x32 shape in {path}; "
            f"found {array.shape}"
        )
    # GLASS convention used by the paper: row 0=ell_eff, row 4=C_ell^BB.
    return {"ell": array[0], "BB": array[4], "source": str(path)}


def _load_pure_covariance(path: str | Path, theta: np.ndarray):
    with np.load(path, allow_pickle=False) as data:
        if "cov_pure_eb" not in data:
            raise KeyError(f"{path} has no cov_pure_eb key; keys={data.files}")
        covariance = _native(data["cov_pure_eb"])
        hartlap_n = hartlap_samples(data)
    if covariance.shape != (120, 120):
        raise ValueError(
            f"Expected 120x120 pure-EB covariance, found {covariance.shape}"
        )
    nbins = len(theta)
    xip = covariance[2 * nbins : 3 * nbins, 2 * nbins : 3 * nbins]
    xim = covariance[3 * nbins : 4 * nbins, 3 * nbins : 4 * nbins]
    cross = covariance[2 * nbins : 3 * nbins, 3 * nbins : 4 * nbins]
    mask = (theta >= PURE_THETA_MIN) & (theta <= PURE_THETA_MAX)
    if int(mask.sum()) != 7:
        warnings.warn(
            f"Pure-EB [12, 83] selection found {mask.sum()} bins, expected 7",
            RuntimeWarning,
        )
    joint = np.block(
        [
            [xip[np.ix_(mask, mask)], cross[np.ix_(mask, mask)]],
            [cross.T[np.ix_(mask, mask)], xim[np.ix_(mask, mask)]],
        ]
    )
    joint_full = np.block([[xip, cross], [cross.T, xim]])
    return {
        "full": covariance,
        "xip_full": xip,
        "xim_full": xim,
        "cross_full": cross,
        "mask": mask,
        "xip_cut": xip[np.ix_(mask, mask)],
        "xim_cut": xim[np.ix_(mask, mask)],
        "cross_cut": cross[np.ix_(mask, mask)],
        "joint_cut": joint,
        "joint_full": joint_full,
        "sigma_xip_full": np.sqrt(np.clip(np.diag(xip), 0, None)),
        "sigma_xim_full": np.sqrt(np.clip(np.diag(xim), 0, None)),
        "hartlap_n": hartlap_n,
    }


def _load_cl_covariance(path: str | Path, ell: np.ndarray, scale=1.0):
    with fits.open(path) as hdul:
        covariance = _native(hdul["COVAR_BB_BB"].data)
    if covariance.shape != (32, 32):
        raise ValueError(f"Expected 32x32 C_l BB covariance, found {covariance.shape}")
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(
            f"C_l covariance scale must be positive and finite, found {scale}"
        )
    # Mock matching assumes the BB covariance is noise dominated, so the
    # stored fiducial BB covariance is scaled by (sigma_e^2 / nbar) ratio^2.
    # Mixed signal-noise pieces would carry one power only; with BB signal=0,
    # the remaining caveat is mask-coupled EE×BB structure.
    covariance *= float(scale)
    mask = ell_bin_mask(ell, ELL_MIN, ELL_MAX)
    if int(mask.sum()) != 17:
        warnings.warn(
            f"Powspace [300, 1600] selection found {mask.sum()} bands, expected 17",
            RuntimeWarning,
        )
    return {
        "full": covariance,
        "scale": float(scale),
        "mask": mask,
        "cut": covariance[np.ix_(mask, mask)],
        "ell_low": get_powspace_bin_edges(ell)[0],
        "ell_high": get_powspace_bin_edges(ell)[1],
    }


def _load_cosebis_covariance(
    path_xi: str | Path,
    path_cov: str | Path,
    *,
    theta_min: float,
    theta_max: float,
):
    """Build the paper-identical COSEBI B covariance on an angular range."""

    gg = treecorr.GGCorrelation(
        min_sep=0.5, max_sep=300.0, nbins=1000, sep_units="arcmin"
    )
    gg.read(str(path_xi))
    start, stop = scale_cut_to_bins(gg, theta_min, theta_max)
    inds = np.arange(start, stop)
    theta_cut = _native(gg.meanr[inds])
    covariance_xipm = _native(np.loadtxt(path_cov))
    nbins = len(gg.meanr)
    if covariance_xipm.shape != (2 * nbins, 2 * nbins):
        raise ValueError(
            f"Expected a {2 * nbins}x{2 * nbins} ξ covariance, "
            f"found {covariance_xipm.shape}"
        )
    covariance_indices = np.concatenate([inds, inds + nbins])
    cosebis = COSEBIS(float(theta_cut.min()), float(theta_cut.max()), 20, precision=120)
    covariance_cosebis = cosebis.cosebis_covariance_from_xipm_covariance(
        theta_cut,
        covariance_xipm[covariance_indices[:, None], covariance_indices],
    )
    covariance_B = covariance_cosebis[20:, 20:]
    sigma = np.sqrt(np.clip(np.diag(covariance_B), 0, None))
    print(
        "COSEBI analytic sigma_n (B, n=1..20): "
        + " ".join(f"{value:.8e}" for value in sigma),
        flush=True,
    )
    return {
        "full": covariance_B,
        "sigma": sigma,
        "theta_cut": theta_cut,
        "requested_cut": (theta_min, theta_max),
    }


def _compute_statistic(data, covariance, *, hartlap_n=None):
    return compute_chi2_pte(data, covariance, n_samples=hartlap_n)


def _statistics_for_mock(cosebis, cosebis_full, pure, cl, covariances):
    bn = cosebis["Bn"]
    bn_full = cosebis_full["Bn"]
    xip = pure["xip_B"]
    xim = pure["xim_B"]
    pure_mask = covariances["pure"]["mask"]
    xip_cut = xip[pure_mask]
    xim_cut = xim[pure_mask]
    xi_total = np.concatenate([xip_cut, xim_cut])
    cl_bb = cl["BB"][covariances["cl"]["mask"]]
    cosebis_cov = covariances["cosebis"]["fiducial"]["full"]
    cosebis_cov_full = covariances["cosebis"]["full_range"]["full"]
    pure_cov = covariances["pure"]
    cl_cov = covariances["cl"]["cut"]

    stat_data = [bn[:6], bn, xip_cut, xim_cut, xi_total, cl_bb]
    stat_cov = [
        cosebis_cov[:6, :6],
        cosebis_cov,
        pure_cov["xip_cut"],
        pure_cov["xim_cut"],
        pure_cov["joint_cut"],
        cl_cov,
    ]
    n_pure = pure_cov["hartlap_n"]
    hartlap = [None, None, n_pure, n_pure, n_pure, None]
    outputs = [
        _compute_statistic(data, covariance, hartlap_n=n_samples)
        for data, covariance, n_samples in zip(stat_data, stat_cov, hartlap)
    ]

    cl_bb_full = cl["BB"]
    full_data = [
        bn_full[:6],
        bn_full,
        xip,
        xim,
        np.concatenate([xip, xim]),
        cl_bb_full,
    ]
    full_cov = [
        cosebis_cov_full[:6, :6],
        cosebis_cov_full,
        pure_cov["xip_full"],
        pure_cov["xim_full"],
        pure_cov["joint_full"],
        covariances["cl"]["full"],
    ]
    full_outputs = [
        _compute_statistic(data, covariance, hartlap_n=n_samples)
        for data, covariance, n_samples in zip(full_data, full_cov, hartlap)
    ]
    return {
        "Bn": bn,
        "Bn_full": bn_full,
        "xip_B_full": xip,
        "xim_B_full": xim,
        "xip_B_cut": xip_cut,
        "xim_B_cut": xim_cut,
        "xi_tot_B": xi_total,
        "cl_BB": cl_bb,
        "cl_BB_full": cl_bb_full,
        "chi2": np.asarray([result[0] for result in outputs], dtype=float),
        "pte": np.asarray([result[1] for result in outputs], dtype=float),
        "dof": np.asarray([result[2] for result in outputs], dtype=int),
        "chi2_full": np.asarray(
            [result[0] for result in full_outputs], dtype=float
        ),
        "pte_full": np.asarray(
            [result[1] for result in full_outputs], dtype=float
        ),
        "dof_full": np.asarray([result[2] for result in full_outputs], dtype=int),
    }


def _json_number(value):
    value = float(value)
    return value if np.isfinite(value) else None


def _json_array(values):
    return [_json_number(value) for value in np.asarray(values).ravel()]


def _sigma_evidence(label, empirical_data, analytic_sigma, labels=None):
    empirical_sigma = np.std(empirical_data, axis=0, ddof=1)
    analytic_sigma = np.asarray(analytic_sigma, dtype=float)
    ratio = empirical_sigma / analytic_sigma
    finite = ratio[np.isfinite(ratio)]
    summary = {
        "min": _json_number(np.min(finite)) if len(finite) else None,
        "max": _json_number(np.max(finite)) if len(finite) else None,
        "median": _json_number(np.median(finite)) if len(finite) else None,
    }
    return {
        "label": label,
        "labels": _json_array(labels) if labels is not None else None,
        "analytic_sigma": _json_array(analytic_sigma),
        "empirical_sigma": _json_array(empirical_sigma),
        "ratio": _json_array(ratio),
        "summary": summary,
    }


def _rank_columns(array):
    return np.apply_along_axis(stats.rankdata, 0, array)


def _correlation_matrix(array, method):
    if array.shape[0] < 2:
        return np.full((array.shape[1], array.shape[1]), np.nan)
    values = _rank_columns(array) if method == "spearman" else array
    return np.corrcoef(values, rowvar=False)


def _chi2_evidence(chi2, pte, dof, pure_hartlap_n):
    result = {}
    for index, name in enumerate(STAT_NAMES):
        values = np.asarray(chi2[:, index], dtype=float)
        ptes = np.asarray(pte[:, index], dtype=float)
        finite = np.isfinite(values) & np.isfinite(ptes)
        values, ptes = values[finite], ptes[finite]
        nu = int(dof[index])
        ks = stats.kstest(ptes, "uniform") if len(ptes) else None
        mean = float(np.mean(values)) if len(values) else np.nan
        variance = float(np.var(values, ddof=1)) if len(values) > 1 else np.nan
        effective = 2.0 * mean**2 / variance if variance > 0 else np.nan
        result[name] = {
            "n": int(len(values)),
            "dof": nu,
            "hartlap_n": pure_hartlap_n if index in (2, 3, 4) else None,
            "chi2_mean": _json_number(mean),
            "chi2_variance": _json_number(variance),
            "expected_mean_nu": nu,
            "expected_variance_2nu": 2 * nu,
            "effective_dof": _json_number(effective),
            "effective_dof_from_mean": _json_number(mean),
            "effective_dof_from_variance": _json_number(variance / 2),
            "ks_pte_uniformity": {
                "statistic": _json_number(ks.statistic) if ks else None,
                "pvalue": _json_number(ks.pvalue) if ks else None,
            },
        }
    return result


def _correlation_evidence(pte, chi2):
    def nested_json(matrix):
        return [[_json_number(value) for value in row] for row in matrix]

    matrices = {}
    for label, values in (("pte", pte), ("chi2", chi2)):
        matrices[label] = {
            "stat_names": list(STAT_NAMES),
            "pearson": _json_array(_correlation_matrix(values, "pearson")),
            "spearman": _json_array(_correlation_matrix(values, "spearman")),
        }
        matrices[label]["shape"] = list(np.asarray(values).shape)
        # The flat arrays above are convenient for strict JSON consumers;
        # retain nested forms for human use and plotting provenance.
        matrices[label]["pearson_matrix"] = nested_json(
            _correlation_matrix(values, "pearson")
        )
        matrices[label]["spearman_matrix"] = nested_json(
            _correlation_matrix(values, "spearman")
        )
    return matrices


def _save_figure(fig, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}", flush=True)


def _plot_sigma_ratios(out_path, stats_payload):
    cosebi_ratio = (
        np.std(stats_payload["Bn"], axis=0, ddof=1) / stats_payload["sigma_cosebis"]
    )
    xip_ratio = (
        np.std(stats_payload["xip_B_full"], axis=0, ddof=1)
        / stats_payload["sigma_xip_full"]
    )
    xim_ratio = (
        np.std(stats_payload["xim_B_full"], axis=0, ddof=1)
        / stats_payload["sigma_xim_full"]
    )
    cl_ratio = (
        np.std(stats_payload["cl_BB"], axis=0, ddof=1) / stats_payload["sigma_cl_cut"]
    )
    theta = stats_payload["theta"]
    ell_cut = stats_payload["ell_cut"]

    # Sampling error of an N-realization standard deviation:
    # sigma(s)/s ~= 1/sqrt(2(N-1)) => each ratio point carries that relative error.
    n_real = stats_payload["Bn"].shape[0]
    rel_err = 1.0 / np.sqrt(2.0 * (n_real - 1))
    eb_kw = dict(marker="o", color=CATEGORICAL_COLORS[0], **SCATTER_ERRORBAR)
    plotted_ratios = (cosebi_ratio, xip_ratio, xim_ratio, cl_ratio)
    plotted_content = np.concatenate(
        [
            np.asarray(ratio)[..., None]
            + np.asarray(ratio * rel_err)[..., None] * np.array([-1.0, 1.0])
            for ratio in plotted_ratios
        ]
    )
    finite_content = plotted_content[np.isfinite(plotted_content)]
    max_dev = np.max(np.abs(finite_content - 1.0))
    shared_ylim = (1.0 - 1.05 * max_dev, 1.0 + 1.05 * max_dev)

    fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.8))
    axes[0].errorbar(
        np.arange(1, 21), cosebi_ratio, yerr=cosebi_ratio * rel_err, **eb_kw
    )
    axes[0].set_xlabel("COSEBI mode $n$")
    axes[0].set_ylabel(r"$\sigma_\mathrm{emp}/\sigma_\mathrm{ana}$")
    axes[0].set_xticks([1, 5, 10, 15, 20])

    pure_theta = np.asarray(theta)
    axes[1].errorbar(
        pure_theta,
        xip_ratio,
        yerr=xip_ratio * rel_err,
        label=r"$\xi_+^B$",
        **eb_kw,
    )
    axes[1].errorbar(
        pure_theta,
        xim_ratio,
        yerr=xim_ratio * rel_err,
        **{**eb_kw, "marker": "s", "color": CATEGORICAL_COLORS[1]},
        label=r"$\xi_-^B$",
    )
    axes[1].axvspan(PURE_THETA_MIN, PURE_THETA_MAX, color="0.85", zorder=0)
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"$\theta$ [arcmin]")
    axes[1].legend(fontsize=6)

    axes[2].errorbar(ell_cut, cl_ratio, yerr=cl_ratio * rel_err, **eb_kw)
    axes[2].axvspan(ELL_MIN, ELL_MAX, color="0.85", zorder=0)
    axes[2].set_xlabel(r"$\ell$")
    axes[2].set_xscale("squareroot" if PAPER_MPLSTYLE else "linear")

    for ax in axes:
        ax.axhline(1.0, color="black", ls="--", lw=0.8)
        ax.grid(alpha=0.2)
        ax.set_ylim(shared_ylim)
    fig.suptitle(
        r"GLASS mock covariance campaign: empirical/analytic $\sigma$", fontsize=9
    )
    _save_figure(fig, out_path)


def _read_data_ptes(path):
    """Read the data's per-statistic PTEs from a min-PTE calibration json.

    Keyed by the same STAT_NAMES the stats payload uses, so the figure can
    look each panel's value up by name.
    """

    payload = json.loads(Path(path).read_text())
    observed = payload["observed_ptes"]
    unknown = set(observed) - set(STAT_NAMES)
    if unknown:
        raise SystemExit(
            "--data-ptes carries statistics this figure does not know: "
            + ", ".join(sorted(unknown))
        )
    return {name: float(value) for name, value in observed.items()}


def _plot_pte_uniformity(out_path, pte, data_ptes=None):
    fig, axes = plt.subplots(2, 3, figsize=(8.0, 5.0), squeeze=False)
    quantiles = (np.arange(1, len(pte) + 1) - 0.5) / len(pte)
    for index, (ax, name) in enumerate(zip(axes.flat, STAT_NAMES)):
        values = np.asarray(pte[:, index])
        ax.hist(values, bins=np.linspace(0, 1, 11), density=True, alpha=0.65)
        ax.set_xlim(0, 1)
        ax.set_xlabel("PTE")
        ax.set_ylabel("density")
        ks = stats.kstest(values, "uniform")
        ax.set_title(f"{STAT_TEX.get(name, name)}\nKS p={ks.pvalue:.3g}", fontsize=8)
        ax.grid(alpha=0.2)
        if data_ptes is not None and name in data_ptes:
            ax.axvline(
                float(data_ptes[name]),
                color="crimson",
                lw=1.2,
                label=f"data: {float(data_ptes[name]):.3g}",
            )
            ax.legend(fontsize=6, loc="best")
        inset = ax.inset_axes([0.66, 0.14, 0.30, 0.24])
        sorted_values = np.sort(values)
        inset.plot(quantiles, sorted_values, ".", ms=1.5)
        inset.plot([0, 1], [0, 1], "k--", lw=0.5)
        inset.set_xlim(0, 1)
        inset.set_ylim(0, 1)
        inset.set_xticks([0, 1])
        inset.set_yticks([0, 1])
        inset.set_title("QQ", fontsize=5)
        inset.tick_params(labelsize=4, length=2, pad=1)
    fig.suptitle("PTE uniformity and uniform QQ diagnostics", fontsize=9)
    fig.tight_layout()
    _save_figure(fig, out_path)


def _heatmap(ax, matrix, title):
    image = ax.imshow(matrix, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(
        range(len(STAT_NAMES)),
        [STAT_TEX.get(n, n) for n in STAT_NAMES],
        rotation=45,
        ha="right",
        fontsize=7,
    )
    ax.set_yticks(
        range(len(STAT_NAMES)), [STAT_TEX.get(n, n) for n in STAT_NAMES], fontsize=7
    )
    for row in range(len(STAT_NAMES)):
        for col in range(len(STAT_NAMES)):
            value = matrix[row, col]
            if np.isfinite(value):
                ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=6)
    ax.set_title(title, fontsize=8)
    return image


def _plot_correlations(out_path, pte, chi2):
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    image = _heatmap(axes[0], _correlation_matrix(pte, "spearman"), "PTE Spearman")
    _heatmap(axes[1], _correlation_matrix(chi2, "spearman"), r"$\chi^2$ Spearman")
    fig.colorbar(image, ax=axes, shrink=0.8, label="correlation")
    fig.suptitle("Cross-statistic correlation of null-test outputs", fontsize=9)
    fig.tight_layout()
    _save_figure(fig, out_path)


def compute_statistics(
    *,
    cosebis_ensemble_path,
    cosebis_ensemble_full_path,
    pure_ensemble_path,
    fine_ensemble_path,
    cl_dir,
    pseudo_cl_dir,
    cosebis_xi_grid,
    cosebis_cov,
    pure_cov,
    cl_cov,
    out_path,
    label="fiducial",
    cl_cov_scale=1.0,
):
    """The compute pass.

    Reads the four mock ensembles and the reference covariance triple, forms
    the six statistics for every realisation at both angular windows, and writes
    one npz holding the per-realisation tables together with everything the
    reporting passes need:
    the analytic σ vectors, the cut grids, and a provenance blob.  Nothing
    downstream of this touches a covariance file.
    """

    cosebis_ensemble = _read_ensemble(cosebis_ensemble_path)
    cosebis_ensemble_full = _read_ensemble(cosebis_ensemble_full_path)
    pure_ensemble = _read_ensemble(pure_ensemble_path)
    fine_ensemble = _read_ensemble(fine_ensemble_path) if fine_ensemble_path else None

    requested_ids = _normalise_ids(cosebis_ensemble["mock_ids"])
    index_of = {
        "cosebis": {mock_id: i for i, mock_id in enumerate(cosebis_ensemble["mock_ids"])},
        "cosebis_full": {
            mock_id: i for i, mock_id in enumerate(cosebis_ensemble_full["mock_ids"])
        },
        "pure": {mock_id: i for i, mock_id in enumerate(pure_ensemble["mock_ids"])},
    }
    if fine_ensemble is not None:
        index_of["fine"] = {
            mock_id: i for i, mock_id in enumerate(fine_ensemble["mock_ids"])
        }

    complete_ids = []
    skipped = {}
    for mock_id in requested_ids:
        missing = []
        if mock_id not in index_of["pure"]:
            missing.append(f"{pure_ensemble['source']}[{mock_id}]")
        if mock_id not in index_of["cosebis_full"]:
            missing.append(f"{cosebis_ensemble_full['source']}[{mock_id}]")
        cl_path = _path_for("cl", cl_dir, mock_id)
        if not cl_path.is_file():
            missing.append(str(cl_path))
        if missing:
            skipped[mock_id] = {"reason": "missing required input", "paths": missing}
        else:
            complete_ids.append(mock_id)
    if not complete_ids:
        raise RuntimeError(
            "No complete mocks found. Required products are the COSEBI and "
            "pure-EB mode ensembles, and the external per-mock cl npy."
        )
    print(
        f"Campaign candidates={len(requested_ids)}, complete={len(complete_ids)}, "
        f"skipped={len(skipped)}",
        flush=True,
    )
    print(
        f"Covariance mode={label}, C_l BB covariance scale={float(cl_cov_scale):.8g}",
        flush=True,
    )

    first_pure = _load_pure(pure_ensemble, index_of["pure"][complete_ids[0]])
    pure_integration_limits = first_pure.get("integration_limits")
    if pure_integration_limits is None:
        warnings.warn(
            "Pure-EB mock vectors have no integration-limits metadata; "
            "their transform window cannot be checked against the stored "
            "paper pure-EB covariance.",
            RuntimeWarning,
        )
    pure_covariances = _load_pure_covariance(pure_cov, first_pure["theta"])
    cosebis_covariances = {
        "fiducial": _load_cosebis_covariance(
            cosebis_xi_grid,
            cosebis_cov,
            theta_min=COSEBI_THETA_MIN,
            theta_max=COSEBI_THETA_MAX,
        ),
        "full_range": _load_cosebis_covariance(
            cosebis_xi_grid,
            cosebis_cov,
            theta_min=FULL_COSEBI_THETA_MIN,
            theta_max=FULL_COSEBI_THETA_MAX,
        ),
    }
    first_cl = _load_cl(_path_for("cl", cl_dir, complete_ids[0]))
    cl_covariances = _load_cl_covariance(cl_cov, first_cl["ell"], scale=cl_cov_scale)
    covariances = {
        "pure": pure_covariances,
        "cosebis": cosebis_covariances,
        "cl": cl_covariances,
    }

    records = []
    stored_scale_metadata = {}
    ell_eff = first_cl["ell"]
    theta = first_pure["theta"]
    for mock_id in complete_ids:
        try:
            cosebis = _load_or_recompute_cosebis(
                cosebis_ensemble,
                fine_ensemble,
                index_of["cosebis"][mock_id],
            )
            cosebis_full = _load_cosebis(
                cosebis_ensemble_full, index_of["cosebis_full"][mock_id]
            )
            pure = _load_pure(pure_ensemble, index_of["pure"][mock_id])
            cl = _load_cl(_path_for("cl", cl_dir, mock_id))
            # theta is the pair-weighted meanr and legitimately varies per mock;
            # the binning (edges) is what must match the reference.
            edge_keys = ("reporting_left_edges", "reporting_right_edges")
            if all(key in first_pure and key in pure for key in edge_keys):
                for edge_key in edge_keys:
                    if not np.allclose(
                        pure[edge_key], first_pure[edge_key], rtol=1e-10
                    ):
                        raise ValueError(
                            f"pure-EB reporting bin edges differ from mock {complete_ids[0]}"
                        )
            elif not np.allclose(pure["theta"], theta, rtol=1e-3):
                raise ValueError(
                    f"pure-EB theta grid differs from mock {complete_ids[0]}"
                )
            if not np.allclose(cl["ell"], ell_eff, rtol=0, atol=1e-8):
                raise ValueError(f"C_l ell grid differs from mock {complete_ids[0]}")
            result = _statistics_for_mock(
                cosebis, cosebis_full, pure, cl, covariances
            )
            if not all(np.all(np.isfinite(result[key])) for key in ("chi2", "pte")):
                raise ValueError("non-finite chi2/PTE")
            records.append(result)
            stored_scale_metadata[mock_id] = {
                "stored_theta_min": cosebis["stored_theta_min"],
                "stored_theta_max": cosebis["stored_theta_max"],
                "action": cosebis["scale_action"],
            }
        except Exception as exc:
            warnings.warn(f"Skipping mock {mock_id}: {exc}", RuntimeWarning)
            skipped[mock_id] = {"reason": str(exc)}

    if not records:
        raise RuntimeError("All candidate mocks failed validation; see warnings above")
    ids = [mock_id for mock_id in complete_ids if mock_id not in skipped]
    arrays = {
        key: np.stack([record[key] for record in records])
        for key in (
            "Bn",
            "Bn_full",
            "xip_B_full",
            "xim_B_full",
            "xip_B_cut",
            "xim_B_cut",
            "xi_tot_B",
            "cl_BB",
        )
    }
    chi2 = np.stack([record["chi2"] for record in records])
    pte = np.stack([record["pte"] for record in records])
    dof = records[0]["dof"]
    pure_mask = pure_covariances["mask"]
    cl_mask = cl_covariances["mask"]

    # The reporting passes are pure functions of what is written here, so the
    # analytic sigma vectors travel with the empirical stack rather than being
    # re-derived from the covariance files a second time.
    provenance = {
        "spec": "mock_covariance_campaign",
        "generated": datetime.now().isoformat(),
        "n_mocks": len(ids),
        "mock_ids": ids,
        "requested_mock_ids": requested_ids,
        "skipped": skipped,
        "ensembles": {
            "cosebis": str(cosebis_ensemble_path),
            "cosebis_full_range": str(cosebis_ensemble_full_path),
            "pure_eb": str(pure_ensemble_path),
            "fine_xi": str(fine_ensemble_path) if fine_ensemble_path else None,
            "cl_dir": str(cl_dir),
        },
        "scale_cuts": {
            "cosebi_requested_exact": [COSEBI_THETA_MIN, COSEBI_THETA_MAX],
            "cosebi_full_range": [FULL_COSEBI_THETA_MIN, FULL_COSEBI_THETA_MAX],
            "pure_xip_and_xim": [PURE_THETA_MIN, PURE_THETA_MAX],
            "pure_full_range_bins": int(len(theta)),
            "cl_bb": [ELL_MIN, ELL_MAX],
            "cl_full_range": [8.0, 2048.0],
            "pure_bins_selected": int(pure_mask.sum()),
            "cl_bands_selected": int(cl_mask.sum()),
        },
        "covariance_mode": {
            "label": str(label),
            "cosebis_covariance": str(cosebis_cov),
            "cosebis_covariance_windows": {
                "fiducial": [COSEBI_THETA_MIN, COSEBI_THETA_MAX],
                "full_range": [FULL_COSEBI_THETA_MIN, FULL_COSEBI_THETA_MAX],
            },
            "pure_eb_covariance": str(pure_cov),
            "pure_eb_hartlap_n": pure_covariances["hartlap_n"],
            "cl_bb_covariance": str(cl_cov),
            "cl_cov_scale": float(cl_cov_scale),
            # The count-ratio constants describe the fallback rescale only, and
            # are omitted when a directly computed mock covariance was used --
            # reporting them there would attach a v1.4.6-suite number to a run
            # that never applied it.
            **(
                {}
                if cl_cov_scale == 1.0
                else {
                    "mockmatching_noise_density_ratio": FIDUCIAL_NEFF / MOCK_NEFF,
                    "mockmatching_noise_variance_ratio": MOCKMATCHED_CL_COV_SCALE,
                    "cl_cov_scale_description": (
                        "Full stored C_l BB covariance multiplied by the supplied "
                        "noise variance ratio. This is exact for the noise-noise "
                        "term only; mixed signal-noise terms would scale by one "
                        "power, but BB has zero signal here, leaving only "
                        "mask-coupled EE×BB caveat terms."
                    ),
                }
            ),
        },
        "cosebi_scale_conventions": stored_scale_metadata,
        "pure_eb_input_grid": {
            "first_mock_integration_limits": pure_integration_limits,
            "covariance_integration_limits": "fine_grid_extent",
            "comparison_warning": pure_integration_limits is None,
        },
    }

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        mock_ids=np.asarray(ids),
        stat_names=np.asarray(STAT_NAMES),
        theta=theta,
        theta_cut=theta[pure_mask],
        ell=ell_eff,
        ell_cut=ell_eff[cl_mask],
        cosebis_theta_cut=np.asarray(
            cosebis_covariances["fiducial"]["theta_cut"]
        ),
        pure_mask=pure_mask,
        cl_mask=cl_mask,
        sigma_cosebis=cosebis_covariances["fiducial"]["sigma"],
        sigma_cosebis_full=cosebis_covariances["full_range"]["sigma"],
        sigma_xip_full=pure_covariances["sigma_xip_full"],
        sigma_xim_full=pure_covariances["sigma_xim_full"],
        sigma_xip_cut=np.sqrt(np.clip(np.diag(pure_covariances["xip_cut"]), 0, None)),
        sigma_xim_cut=np.sqrt(np.clip(np.diag(pure_covariances["xim_cut"]), 0, None)),
        sigma_cl_cut=np.sqrt(np.clip(np.diag(cl_covariances["cut"]), 0, None)),
        sigma_cl_full=np.sqrt(np.clip(np.diag(cl_covariances["full"]), 0, None)),
        **arrays,
        chi2=chi2,
        pte=pte,
        dof=dof,
        chi2_full=np.stack([record["chi2_full"] for record in records]),
        pte_full=np.stack([record["pte_full"] for record in records]),
        dof_full=records[0]["dof_full"],
        provenance=np.asarray(json.dumps(provenance)),
    )
    print(f"Saved {out_path}", flush=True)

    # One-mock sanity gate: it is diagnostic rather than a campaign veto,
    # since a single Gaussian realization can naturally land below PTE=0.05.
    first = records[0]
    finite = all(
        np.all(np.isfinite(first[key]))
        for key in ("chi2", "pte", "xip_B_cut", "xim_B_cut")
    )
    xip_zmax = float(
        np.max(
            np.abs(first["xip_B_cut"]) / np.sqrt(np.diag(pure_covariances["xip_cut"]))
        )
    )
    xim_zmax = float(
        np.max(
            np.abs(first["xim_B_cut"]) / np.sqrt(np.diag(pure_covariances["xim_cut"]))
        )
    )
    sanity_pass = finite and xip_zmax <= 2.5 and xim_zmax <= 2.5
    print(
        "One-mock sanity: "
        f"finite={finite}, max |xi+^B|/sigma={xip_zmax:.2f}, "
        f"max |xi-^B|/sigma={xim_zmax:.2f}, pass(~2sigma)={sanity_pass}",
        flush=True,
    )
    if not sanity_pass:
        warnings.warn(
            "One-mock pure-B sanity gate is outside the loose 2.5σ diagnostic",
            RuntimeWarning,
        )

    # Provenance check: compare the first
    # dozen IDs for which the in-repo FITS product exists; the external array's
    # row 4 is the BB vector.
    compared = 0
    for mock_id in ids:
        if compared >= 12:
            break
        pseudo_path = _path_for("pseudo_cl", pseudo_cl_dir, mock_id)
        if not pseudo_path.is_file():
            continue
        with fits.open(pseudo_path) as hdul:
            pseudo_bb = _native(hdul["PSEUDO_CELL"].data["BB"])
        external_bb = _native(
            np.load(_path_for("cl", cl_dir, mock_id), allow_pickle=False)[4]
        )
        if pseudo_bb.shape != external_bb.shape:
            print(
                f"WARNING C_l provenance {mock_id}: shape {external_bb.shape} vs {pseudo_bb.shape}",
                flush=True,
            )
            continue
        difference = external_bb - pseudo_bb
        scale = np.maximum(np.abs(pseudo_bb), 1e-30)
        max_abs = float(np.max(np.abs(difference)))
        max_relative = float(np.max(np.abs(difference) / scale))
        rms = float(np.sqrt(np.mean(difference**2)))
        agreement = max_relative <= 1e-3 and np.allclose(
            external_bb, pseudo_bb, rtol=1e-5, atol=1e-14
        )
        print(
            f"C_l provenance {mock_id}: max_abs={max_abs:.4e}, rms={rms:.4e}, "
            f"max_relative={max_relative:.4e}, allclose={agreement}",
            flush=True,
        )
        if max_relative > 1e-3:
            print(
                f"WARNING: material external-vs-FITS BB disagreement for mock {mock_id}",
                flush=True,
            )
        compared += 1

    return out_path


# ── The reporting passes ────────────────────────────────────────────────────


def read_statistics(path):
    """Load a compute-pass npz, with its provenance blob decoded."""

    with np.load(path, allow_pickle=False) as data:
        payload = {key: data[key] for key in data.files}
    payload["provenance"] = json.loads(str(payload["provenance"]))
    payload["stat_names"] = [str(name) for name in payload["stat_names"]]
    return payload


def emit_evidence(payload, out_path):
    """The campaign's machine-readable result, from the stored stack alone."""

    sigma_evidence = {
        "COSEBI_B_modes_1_to_20": _sigma_evidence(
            "COSEBI B modes 1--20",
            payload["Bn"],
            payload["sigma_cosebis"],
            labels=np.arange(1, 21),
        ),
        "xi_plus_B_20_bins": _sigma_evidence(
            "ξ+^B, all 20 reporting bins",
            payload["xip_B_full"],
            payload["sigma_xip_full"],
            labels=payload["theta"],
        ),
        "xi_minus_B_20_bins": _sigma_evidence(
            "ξ−^B, all 20 reporting bins",
            payload["xim_B_full"],
            payload["sigma_xim_full"],
            labels=payload["theta"],
        ),
        "xi_plus_B_7_cut_bins": _sigma_evidence(
            "ξ+^B, 12--83 arcmin cut",
            payload["xip_B_cut"],
            payload["sigma_xip_cut"],
            labels=payload["theta_cut"],
        ),
        "xi_minus_B_7_cut_bins": _sigma_evidence(
            "ξ−^B, 12--83 arcmin cut",
            payload["xim_B_cut"],
            payload["sigma_xim_cut"],
            labels=payload["theta_cut"],
        ),
        "C_ell_BB_17_cut_bands": _sigma_evidence(
            "C_l^BB, 17 bands fully inside 300--1600",
            payload["cl_BB"],
            payload["sigma_cl_cut"],
            labels=payload["ell_cut"],
        ),
    }
    evidence = {
        **payload["provenance"],
        "sigma_ratios": sigma_evidence,
        "chi2_calibration": _chi2_evidence(
            payload["chi2"],
            payload["pte"],
            payload["dof"],
            payload["provenance"]["covariance_mode"].get(
                "pure_eb_hartlap_n", MC_SAMPLES_WITHOUT_RECORD
            ),
        ),
        "correlations": _correlation_evidence(payload["pte"], payload["chi2"]),
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as stream:
        json.dump(evidence, stream, indent=2, allow_nan=False)
    print(f"Saved {out_path}", flush=True)
    return out_path


EMITTERS = {
    "evidence": emit_evidence,
    "sigma-ratio-figure": lambda payload, out_path: _plot_sigma_ratios(
        out_path, payload
    ),
    "pte-uniformity-figure": lambda payload, out_path, data_ptes=None: (
        _plot_pte_uniformity(out_path, payload["pte"], data_ptes)
    ),
    "correlation-figure": lambda payload, out_path: _plot_correlations(
        out_path, payload["pte"], payload["chi2"]
    ),
}


def _parser(argv=None):
    defaults = _default_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--emit-stats",
        help="Compute pass: write the per-realisation statistics npz here.",
    )
    parser.add_argument(
        "--from-stats",
        help="Reporting pass: read a compute-pass npz instead of recomputing.",
    )
    parser.add_argument(
        "--emit",
        choices=tuple(EMITTERS),
        help="Which derived product the reporting pass writes.",
    )
    parser.add_argument(
        "--out-file",
        help="Exact path the selected --emit product is written to.",
    )
    parser.add_argument(
        "--data-ptes",
        help=(
            "Optional mock_min_pte_calibration.json.  Its `observed_ptes` — the "
            "real catalogue's PTE per statistic — are drawn as a red vertical "
            "line on the pte-uniformity figure.  Omit it and the figure is "
            "drawn from the mocks alone."
        ),
    )
    parser.add_argument("--cosebis-ensemble", help="Fiducial-window COSEBI mode ensemble npz.")
    parser.add_argument(
        "--cosebis-ensemble-full",
        help="Full-range COSEBI mode ensemble npz.",
    )
    parser.add_argument("--pure-ensemble", help="Pure-mode ensemble npz.")
    parser.add_argument(
        "--fine-ensemble",
        help=(
            "Fine ξ± ensemble npz, used only to recompute COSEBIs whose stored "
            "θ cut is not the fiducial one."
        ),
    )
    parser.add_argument("--cl-dir", default=defaults["cl_dir"])
    parser.add_argument("--pseudo-cl-dir", default=defaults["pseudo_cl_dir"])
    parser.add_argument("--cosebis-xi-grid", default=defaults["cosebis_xi_grid"])
    parser.add_argument("--cosebis-cov", default=defaults["cosebis_cov"])
    parser.add_argument("--pure-cov", default=defaults["pure_cov"])
    parser.add_argument("--cl-cov", default=defaults["cl_cov"])
    parser.add_argument(
        "--cl-cov-scale",
        type=float,
        default=1.0,
        help="Multiplicative scale for the stored C_l BB covariance (default: 1).",
    )
    parser.add_argument(
        "--label",
        default="fiducial",
        help="Covariance suite label written to evidence.json (e.g. fiducial/mockmatched).",
    )
    parser.add_argument(
        "--reference-covariance",
        choices=("data_covariance", "mockmatched"),
        default="data_covariance",
        help=(
            "Which analytic covariance the mock scatter is compared against: "
            "'data_covariance' = the paper's own covariance at the real catalog's "
            "n_eff (--cosebis-cov/--pure-cov/--cl-cov as given); "
            "'mockmatched' = the same CosmoCov/MC chain regenerated at the "
            "mocks' n_eff (--cosebis-cov-mockmatched/--pure-cov-mockmatched), "
            "with the stored C_l^BB covariance rescaled by the noise-variance "
            f"ratio {MOCKMATCHED_CL_COV_SCALE:.6g}."
        ),
    )
    parser.add_argument(
        "--scale-cut-convention",
        choices=tuple(COSEBI_THETA_BY_CONVENTION),
        default="grid_bin_edges",
        help=(
            "How the fiducial COSEBI θ cut is read: 'grid_bin_edges' (what the "
            "pipeline cut at) or 'rounded_arcmin' (what the paper text prints)."
        ),
    )
    parser.add_argument(
        "--cl-cov-mockmatched",
        help=(
            "C_l^BB covariance computed directly on a mock catalogue "
            "(mockmatched only).  When given it is used unscaled; without it "
            "the data covariance is rescaled by the squared count ratio."
        ),
    )
    parser.add_argument(
        "--cosebis-cov-mockmatched",
        help="ξ± covariance regenerated at the mocks' galaxy count (mockmatched only).",
    )
    parser.add_argument(
        "--pure-cov-mockmatched",
        help="Pure-EB MC covariance regenerated at the mocks' galaxy count.",
    )
    return parser.parse_args(argv)


def _resolve_reference_covariance(args):
    """Return (cosebis_cov, pure_cov, cl_cov, cl_cov_scale, label) for the selection."""

    if args.reference_covariance == "data_covariance":
        return (
            args.cosebis_cov,
            args.pure_cov,
            args.cl_cov,
            args.cl_cov_scale,
            args.label,
        )
    missing = [
        name
        for name, value in (
            ("--cosebis-cov-mockmatched", args.cosebis_cov_mockmatched),
            ("--pure-cov-mockmatched", args.pure_cov_mockmatched),
        )
        if not value
    ]
    if missing:
        raise SystemExit(
            "--reference-covariance mockmatched requires " + ", ".join(missing)
        )
    # A directly computed mock C_l^BB covariance, where one exists, replaces the
    # count-ratio rescale of the data covariance: it is a measurement of the
    # mock's own noise rather than an extrapolation from the data's.
    if args.cl_cov_mockmatched:
        cl_cov, cl_cov_scale = args.cl_cov_mockmatched, 1.0
    else:
        cl_cov, cl_cov_scale = args.cl_cov, MOCKMATCHED_CL_COV_SCALE
    return (
        args.cosebis_cov_mockmatched,
        args.pure_cov_mockmatched,
        cl_cov,
        cl_cov_scale,
        args.label,
    )


def _apply_scale_cut_convention(convention):
    """Rebind the θ-cut globals to the selected convention.

    The script reads COSEBI_THETA_MIN/MAX as module constants in six
    places and PURE_THETA_MIN/MAX in three; rebinding them once here keeps the
    decision explicit at the CLI without threading a parameter through the
    whole gather.  The pure-mode cut selects the same seven reporting bins
    either way (see PURE_THETA_BY_CONVENTION), so this moves only the boundary
    values recorded in the evidence file and shaded in the figures.
    """

    global COSEBI_THETA_MIN, COSEBI_THETA_MAX, PURE_THETA_MIN, PURE_THETA_MAX
    COSEBI_THETA_MIN, COSEBI_THETA_MAX = COSEBI_THETA_BY_CONVENTION[convention]
    PURE_THETA_MIN, PURE_THETA_MAX = PURE_THETA_BY_CONVENTION[convention]


def _standalone_main(argv=None):
    args = _parser(argv)
    _apply_scale_cut_convention(args.scale_cut_convention)

    if args.emit_stats:
        missing = [
            name
            for name, value in (
                ("--cosebis-ensemble", args.cosebis_ensemble),
                ("--cosebis-ensemble-full", args.cosebis_ensemble_full),
                ("--pure-ensemble", args.pure_ensemble),
            )
            if not value
        ]
        if missing:
            raise SystemExit("--emit-stats requires " + ", ".join(missing))
        cosebis_cov, pure_cov, cl_cov, cl_cov_scale, label = (
            _resolve_reference_covariance(args)
        )
        compute_statistics(
            cosebis_ensemble_path=args.cosebis_ensemble,
            cosebis_ensemble_full_path=args.cosebis_ensemble_full,
            pure_ensemble_path=args.pure_ensemble,
            fine_ensemble_path=args.fine_ensemble,
            cl_dir=args.cl_dir,
            pseudo_cl_dir=args.pseudo_cl_dir,
            cosebis_xi_grid=args.cosebis_xi_grid,
            cosebis_cov=cosebis_cov,
            pure_cov=pure_cov,
            cl_cov=cl_cov,
            out_path=args.emit_stats,
            label=label,
            cl_cov_scale=cl_cov_scale,
        )
        return

    if not (args.from_stats and args.emit and args.out_file):
        raise SystemExit(
            "give either --emit-stats, or --from-stats with --emit and --out-file"
        )
    extra = {}
    if args.emit == "pte-uniformity-figure" and args.data_ptes:
        extra["data_ptes"] = _read_data_ptes(args.data_ptes)
    EMITTERS[args.emit](read_statistics(args.from_stats), args.out_file, **extra)


if __name__ == "__main__":
    _standalone_main()
