"""The COSEBI transform of ξ_± and of its covariance, for one angular range.

Both act on the same bins of the integration grid and the same filters, built
over the bin centres inside the cut: E_n and B_n from the measured ξ_±, and
Tᵀ C_ξ T from the ξ_± covariance and those bin centres alone.
"""

import numpy as np
from cosmo_numba.B_modes.cosebis import COSEBIS

from sp_validation.b_modes import scale_cut_to_bins


def _filters(gg, nmodes, scale_cut):
    """The COSEBI filters of a cut and the integration-grid bins inside it."""
    start_bin, stop_bin = scale_cut_to_bins(gg, *scale_cut)
    inds = np.arange(start_bin, stop_bin)
    theta_cut = gg.meanr[inds]
    # precision=120 avoids sympy root convergence failures for modes 11 and up.
    cosebis = COSEBIS(
        theta_min=np.min(theta_cut),
        theta_max=np.max(theta_cut),
        N_max=nmodes,
        precision=120,
    )
    return cosebis, inds


def cosebis_modes(gg, nmodes, scale_cut):
    """E_n and B_n of a TreeCorr GGCorrelation over one cut."""
    cosebis, inds = _filters(gg, nmodes, scale_cut)
    return cosebis.cosebis_from_xipm(gg.meanr[inds], gg.xip[inds], gg.xim[inds], parallel=True)


def cosebis_covariance(gg, cov_xipm, nmodes, scale_cut):
    """The (2 nmodes)² E_n ⊕ B_n covariance over one cut, from the ξ_± covariance."""
    cosebis, inds = _filters(gg, nmodes, scale_cut)
    cov_inds = np.concatenate([inds, inds + len(gg.meanr)])
    return cosebis.cosebis_covariance_from_xipm_covariance(
        gg.meanr[inds], cov_xipm[cov_inds[:, None], cov_inds]
    )
