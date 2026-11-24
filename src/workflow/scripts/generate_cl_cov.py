"""Gaussian (iNKA) covariance of the catalog-based Cℓ bandpowers.

The covariance of one catalog version, from its entry in the catalog summary
(``catalog_summary.py``) rather than the catalog itself: the weighted
nside-1024 count map is the window, and the analytic shot-noise bias is
Ω_pix · Σ w²(e1² + e2²)/2 / N_pix, the mean of the per-pixel variance map. The
theory spectrum is CCL's Cℓ for the version's n(z) at the requested blind,
multiplied by the squared pixel window, plus that noise, coupled through the
mode-coupling matrix. This is the computation of sp_validation's
``CosmologyValidation.calculate_pseudo_cl_eb_cov`` with the analytic noise
bias and coupled iNKA input, which this script reproduces.

    python generate_cl_cov.py --ver SP_v1.4.6.3_leak_corr --summary summary.npz \
        --cat-config cat_config.yaml --nbins 32 --blind blind_a --out-file cov.fits

``--ver`` names a cat_config entry, with ``_leak_corr`` for the
leakage-corrected ellipticities; its summary key is the entry without that
suffix.
"""

import argparse
import json
import re
import tempfile

import healpy as hp
import numpy as np
import pymaster as nmt
import yaml
from astropy.io import fits
from sp_validation.cosmo_val import CosmologyValidation
from sp_validation.cosmology import get_theo_c_ell

from catalog_summary import load, window
from decision_values import FIDUCIAL_BLIND

LMIN = 8


def cl_covariance(version, summary, cat_config, nside, blind, cosmo_params, binning, nbins, power):
    """The 16 EE/EB/BE/BB covariance blocks, as an HDUList."""
    base = version.removesuffix("_leak_corr")
    with open(cat_config) as f:
        cc = yaml.safe_load(f)
    nz_path = cc[base]["shear"]["redshift_path"]
    if blind:
        nz_path = re.sub(r"_[ABC]\.txt$", f"_{blind}.txt", nz_path)
    if cosmo_params:
        cosmo_params = {
            "Omega_m": cosmo_params["Omega_m"],
            "Omega_b": cosmo_params["Omega_b"],
            "h": cosmo_params["h"],
            "sig8": cosmo_params["sigma_8"],
            "ns": cosmo_params["n_s"],
        }
    # Supplies the cosmology and the NaMaster binning and field helpers.
    with tempfile.TemporaryDirectory(prefix="cl_cov_") as scratch:
        cv = CosmologyValidation(
            versions=[version],
            catalog_config=cat_config,
            output_dir=scratch,
            binning=binning,
            nside=nside,
            cell_method="catalog",
            cosmo_params=cosmo_params,
            n_ell_bins=nbins,
            power=power,
        )

    entry = load(summary, base)
    n_gal = window(entry)
    s = entry["sum_w2e2_leak_corr" if version.endswith("_leak_corr") else "sum_w2e2"]
    noise_bias = hp.nside2pixarea(nside) * s / hp.nside2npix(nside)

    lmax = 2 * nside
    ell = np.arange(1, lmax + 1)
    pw = hp.pixwin(nside, lmax=lmax)[1 : lmax + 1]
    z, dndz = np.loadtxt(nz_path, unpack=True)
    print(f"Theory n(z): {nz_path}")
    fiducial_cl = get_theo_c_ell(ell=ell, z=z, nz=dndz, backend="ccl", cosmo=cv.cosmo) * pw**2

    b = cv.get_namaster_bin(LMIN, lmax, lmax - 1)
    f, wsp = cv.get_field_and_workspace_from_map(n_gal, lmax - 1, b)

    noise_bias_cl = np.zeros((4, lmax))
    noise_bias_cl[0, :] = noise_bias
    noise_bias_cl[3, :] = noise_bias
    noise_bias_cl = b.unbin_cell(wsp.decouple_cell(noise_bias_cl))
    lowest_ell = b.get_ell_list(0)[0]
    noise_bias_cl[:, :lowest_ell] = noise_bias_cl[:, [lowest_ell]]

    zero = 0.0 * fiducial_cl
    fiducial_cl = np.array([fiducial_cl, zero, zero, zero]) + noise_bias_cl
    coupling = np.reshape(wsp.get_coupling_matrix(), (4, lmax, 4, lmax), order="F")
    fiducial_cl = np.tensordot(coupling, fiducial_cl) / np.mean(n_gal**2)

    cw = nmt.NmtCovarianceWorkspace.from_fields(f, f, f, f)
    n_ell = b.get_n_bands()
    cov = nmt.gaussian_covariance(
        cw, 2, 2, 2, 2, fiducial_cl, fiducial_cl, fiducial_cl, fiducial_cl, wsp, wb=wsp
    ).reshape([n_ell, 4, n_ell, 4])

    pols = ["EE", "EB", "BE", "BB"]
    hdu = fits.HDUList()
    for i, pa in enumerate(pols):
        for j, pb in enumerate(pols):
            hdu.append(fits.ImageHDU(cov[:, i, :, j], name=f"COVAR_{pa}_{pb}"))
    return hdu


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ver", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--cat-config", required=True)
    ap.add_argument("--out-file", required=True)
    ap.add_argument("--nside", type=int, default=1024)
    ap.add_argument("--binning", choices=["linear", "logspace", "powspace"], default="powspace")
    ap.add_argument("--nbins", type=int, required=True)
    ap.add_argument("--power", type=float, default=0.5)
    ap.add_argument("--blind", type=FIDUCIAL_BLIND.__getitem__)
    ap.add_argument("--cosmo-json")
    a = ap.parse_args(argv)
    cosmo_params = None
    if a.cosmo_json:
        with open(a.cosmo_json) as f:
            cosmo_params = json.load(f)
    hdu = cl_covariance(a.ver, a.summary, a.cat_config, a.nside, a.blind, cosmo_params, a.binning, a.nbins, a.power)
    hdu.writeto(a.out_file, overwrite=True)
    print(f"Saved to: {a.out_file}")


if __name__ == "__main__":
    main()
