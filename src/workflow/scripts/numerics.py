"""Null-test statistics and bandpower binning shared by the analysis scripts.

The chi-squared and PTE of a null test (with the Hartlap factor), the NaMaster
square-root bandpower edges and the ell windows they select, and the PTE
threshold each decision option stands for.
"""

import numpy as np
from scipy import stats


PTE_THRESHOLD_VALUES = {
    "nominal_0p05": 0.05,
    "calibrated_0p01": 0.01,
}


def resolve_pte_threshold(value):
    """Resolve an ASTRA PTE-threshold option ID or numeric value."""
    if value in PTE_THRESHOLD_VALUES:
        threshold = PTE_THRESHOLD_VALUES[value]
    else:
        threshold = float(value)
    if not 0 < threshold < 0.5:
        raise ValueError(f"PTE threshold must be between 0 and 0.5, got {threshold}")
    return threshold


def compute_chi2_pte(data, covariance, n_samples=None):
    """Compute chi-squared and PTE for null test.

    Parameters
    ----------
    data : array_like
        Data vector (e.g., B-mode signal).
    covariance : array_like
        Covariance matrix.
    n_samples : int, optional
        Number of samples used to estimate the covariance (e.g., MC samples
        or jackknife patches). When provided, the Hartlap correction factor
        (n_samples - n_bins - 2) / (n_samples - 1) is applied to debias the
        inverse covariance estimate. Leave as None for analytical covariances.

    Returns
    -------
    chi2 : float
        Chi-squared value.
    pte : float
        Probability to exceed (survival function).
    dof : int
        Degrees of freedom (length of data).
    """
    chi2 = float(data @ np.linalg.solve(covariance, data))
    dof = len(data)
    if n_samples is not None:
        hartlap_factor = (n_samples - dof - 2) / (n_samples - 1)
        chi2 *= hartlap_factor
    pte = stats.chi2.sf(chi2, dof)
    return chi2, pte, dof


def get_powspace_bin_edges(ell_eff, lmin=8, lmax=2048, power=0.5):
    """Reconstruct NaMaster powspace bin edges from effective ells.

    Reproduces the binning logic in CosmologyValidation.get_namaster_bin():
    edges are placed at ell^power linspace, then each multipole is assigned
    to a bin via np.digitize.

    Parameters
    ----------
    ell_eff : array_like
        Effective (mean) ell per bin, as stored in Cl FITS files.
    lmin, lmax : int
        Multipole range used in the NaMaster binning.
    power : float
        Exponent for power-law spacing (0.5 = sqrt).

    Returns
    -------
    ell_low : ndarray
        Lower multipole edge (inclusive) for each bin.
    ell_high : ndarray
        Upper multipole edge (inclusive) for each bin.
    """
    n_ell_bins = len(ell_eff)
    ells = np.arange(lmin, lmax + 1)
    start = np.power(lmin, power)
    end = np.power(lmax, power)
    bins_ell = np.power(np.linspace(start, end, n_ell_bins + 1), 1 / power)
    bpws = np.digitize(ells.astype(float), bins_ell) - 1
    bpws[0] = 0
    bpws[-1] = n_ell_bins - 1

    ell_low = np.array([ells[bpws == i][0] for i in range(n_ell_bins)])
    ell_high = np.array([ells[bpws == i][-1] for i in range(n_ell_bins)])
    return ell_low, ell_high


def ell_bin_mask(ell_eff, ell_min, ell_max, **kwargs):
    """Return boolean mask selecting bins fully within [ell_min, ell_max].

    A bin is included only if its lower edge >= ell_min AND its upper
    edge <= ell_max.

    Parameters
    ----------
    ell_eff : array_like
        Effective ell per bin.
    ell_min, ell_max : float
        Scale cut boundaries.
    **kwargs
        Passed to get_powspace_bin_edges.

    Returns
    -------
    mask : ndarray of bool
    """
    ell_low, ell_high = get_powspace_bin_edges(ell_eff, **kwargs)
    return (ell_low >= ell_min) & (ell_high <= ell_max)


def ell_bin_index(ell_eff, ell_min, ell_max, **kwargs):
    """Return (i_min, i_max) indices of first/last bins fully within cut.

    Semantics match the PTE matrix convention: pte_matrix[i_min, i_max]
    uses slice(i_min, i_max + 1), so i_max is inclusive.

    Parameters
    ----------
    ell_eff : array_like
        Effective ell per bin.
    ell_min, ell_max : float
        Scale cut boundaries.
    **kwargs
        Passed to get_powspace_bin_edges.

    Returns
    -------
    i_min, i_max : int
    """
    mask = ell_bin_mask(ell_eff, ell_min, ell_max, **kwargs)
    indices = np.where(mask)[0]
    return int(indices[0]), int(indices[-1])
