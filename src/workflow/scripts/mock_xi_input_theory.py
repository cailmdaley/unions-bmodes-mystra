"""The theory xi_pm the GLASS suite was generated from, on the reporting grid.

The mocks are lognormal realisations built by Sacha Guerrini's
``make_unions_glass_sim.py``: CAMB matter shells under Planck-18 with HMcode
``mead2020_feedback`` (log T_AGN = 7.8) and As rescaled to sigma_8 = 0.8102,
lensed through ``glass.MultiPlaneConvergence`` onto an NSide-4096 map.  No
theory xi_pm ships with the suite, so this script rebuilds the generator's
cosmology verbatim from the project's own ``cosmology_snapshot`` -- the same
Planck-18 numbers the generator reads out of ``astropy`` -- takes the lensing
source C_ell for the fiducial n(z) from CAMB, and transforms it to xi_pm.

The transform is truncated at the simulation's own ``lmax = nside``: the maps
carry no power above it, so an untruncated theory would sit above the mock mean
at the smallest separations for a reason that has nothing to do with the mocks.
Both are written -- ``xip``/``xim`` truncated, ``xip_untruncated``/
``xim_untruncated`` not -- so the size of that effect is readable rather than
assumed.
"""

from __future__ import annotations

import argparse
import json
import tempfile

import camb
import numpy as np
import pyccl as ccl


def build_camb_params(snapshot, log_t_agn):
    """Replicate ``make_unions_glass_sim.Sky.__init__``'s CAMB setup."""
    h = snapshot["h"]
    omega_m = snapshot["Omega_m"]
    omega_b = snapshot["Omega_b"]
    n_s = snapshot["n_s"]
    m_nu = snapshot["m_nu"]

    pars = camb.CAMBparams()
    # Two passes, as the generator does: the first only to read omeganu back
    # out, the second to take it out of the cold-dark-matter density.
    pars.set_cosmology(H0=100 * h, omch2=(omega_m - omega_b) * h**2,
                       ombh2=omega_b * h**2, mnu=m_nu)
    omega_c = omega_m - omega_b - pars.omeganu
    pars.set_cosmology(H0=100 * h, omch2=omega_c * h**2,
                       ombh2=omega_b * h**2, mnu=m_nu)

    a_s = snapshot["A_s"]
    pars.InitPower.set_params(As=a_s, ns=n_s)
    pars.WantTransfer = True
    pars.NonLinear = camb.model.NonLinear_both
    pars.NonLinearModel.set_params(halofit_version="mead2020_feedback",
                                   HMCode_logT_AGN=log_t_agn)
    pars.set_matter_power(nonlinear=True, kmax=20)

    # As is a placeholder until sigma_8 is imposed, exactly as in the generator.
    sigma8_temp = float(camb.get_results(pars).get_sigma8_0())
    a_s *= (snapshot["sigma_8"] / sigma8_temp) ** 2
    init_power = camb.InitialPowerLaw()
    init_power.set_params(As=a_s, ns=n_s)
    pars.InitPower = init_power
    sigma8 = float(camb.get_results(pars).get_sigma8_0())
    if not np.isclose(sigma8, snapshot["sigma_8"], rtol=1e-4):
        raise SystemExit(f"As rescaling failed: sigma_8 = {sigma8}")
    return pars, a_s, sigma8


def load_nz(version, blind, cat_config):
    """The (version, blind) n(z), through sp_validation's own resolver.

    Same call ``mock_matched_pure_eb_cov.py`` makes, so the redshift
    distribution here is the one every other stage of the analysis uses.
    """
    from sp_validation.cosmo_val import CosmologyValidation

    validation = CosmologyValidation(
        versions=[version],
        catalog_config=cat_config,
        output_dir=tempfile.mkdtemp(prefix="mock_xi_input_theory_"),
    )
    validation.blind = blind
    z, nz = validation.get_redshift(version)
    return np.asarray(z, dtype=float), np.asarray(nz, dtype=float)


def source_cls(pars, z, nz, lmax, zmax):
    """C_ell^kappa for the lensing source window, CAMB's raw convention."""
    mask = z <= zmax
    pars.SourceWindows = [
        camb.sources.SplinedSourceWindow(z=z[mask], W=nz[mask],
                                         source_type="lensing")
    ]
    pars.set_for_lmax(lmax)
    pars.SourceTerms.limber_windows = True
    results = camb.get_results(pars)
    cls = results.get_source_cls_dict(lmax=lmax, raw_cl=True)
    return np.asarray(cls["W1xW1"], dtype=float)


def bessel_sum(ell, cl, theta_arcmin):
    """xi_pm as an explicit sum over the C_ell that are actually supplied.

    xi_pm(theta) = (1 / 2 pi) sum_l l C_l J_{0,4}(l theta).  The flat-sky
    kernel costs ~0.1% against the curved-sky Legendre sum on these
    separations, which is the price of the property that matters here: the
    band limit is *honest*.  CCL's FFTLog path extrapolates the input
    spectrum in log-log beyond its last multipole, so a "truncated" transform
    run through it silently reinstates the power the simulation does not
    have -- and xi_- at 1 arcmin is built almost entirely out of l ~ 10^4.
    """
    from scipy.special import jv

    ell = np.asarray(ell, dtype=float)
    cl = np.asarray(cl, dtype=float)
    theta = np.asarray(theta_arcmin, dtype=float) * (np.pi / 180.0 / 60.0)
    x = np.outer(theta, ell)
    weight = ell * cl / (2.0 * np.pi)
    return (jv(0, x) @ weight, jv(4, x) @ weight)


def curved_sky(ell, cl, theta_arcmin):
    """The unbanded curved-sky transform, for the flat-sky cross-check."""
    theta_deg = np.asarray(theta_arcmin, dtype=float) / 60.0
    cosmo = ccl.CosmologyVanillaLCDM()
    kw = dict(cosmo=cosmo, ell=ell, C_ell=cl, theta=theta_deg,
              method="fftlog")
    return (ccl.correlation(type="GG+", **kw),
            ccl.correlation(type="GG-", **kw))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cosmology-json", required=True)
    parser.add_argument("--reporting-ensemble", required=True)
    parser.add_argument("--cat-config", required=True)
    parser.add_argument("--version", default="SP_v1.4.6.3_leak_corr")
    parser.add_argument("--blind", default="A")
    parser.add_argument("--nside", type=int, default=4096,
                        help="The suite's HEALPix NSide; lmax = nside.")
    parser.add_argument("--zmax", type=float, default=3.0,
                        help="The generator's redshift cut on the n(z).")
    parser.add_argument("--lmax-untruncated", type=int, default=20000,
                        help="Wide multipole ceiling the band limit is "
                             "measured against.")
    parser.add_argument("--log-t-agn", type=float, default=7.8)
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args()

    snapshot = json.load(open(args.cosmology_json))
    ensemble = np.load(args.reporting_ensemble, allow_pickle=False)
    # The bin centres the mocks actually report at: TreeCorr's pair-weighted
    # mean separation, averaged over the suite.  Evaluating theory at the
    # nominal geometric centre instead would show up as a percent-level tilt.
    theta = np.asarray(ensemble["meanr"], dtype=float).mean(axis=0)

    pars, a_s, sigma8 = build_camb_params(snapshot, args.log_t_agn)
    z, nz = load_nz(args.version, args.blind, args.cat_config)
    lmax = args.nside

    # One spectrum, to the wide lmax; the band limit is then applied to the
    # transform, so truncated and untruncated differ in nothing else.
    cl_long = source_cls(pars, z, nz, args.lmax_untruncated, args.zmax)
    ell_long = np.arange(cl_long.size, dtype=float)
    band = (ell_long >= 2) & (ell_long <= lmax)
    wide = ell_long >= 2

    xip, xim = bessel_sum(ell_long[band], cl_long[band], theta)
    xip_u, xim_u = bessel_sum(ell_long[wide], cl_long[wide], theta)
    xip_curved, xim_curved = curved_sky(ell_long[wide], cl_long[wide], theta)
    cl = cl_long[: lmax + 1]

    np.savez(
        args.out_file,
        theta=theta,
        xip=xip,
        xim=xim,
        xip_untruncated=xip_u,
        xim_untruncated=xim_u,
        xip_curved_sky=xip_curved,
        xim_curved_sky=xim_curved,
        ell=np.arange(cl.size, dtype=float),
        cl_kappa=cl,
        z=z,
        nz=nz,
        provenance=json.dumps({
            "cosmology": snapshot,
            "A_s_rescaled": a_s,
            "sigma_8_realised": sigma8,
            "halofit_version": "mead2020_feedback",
            "HMCode_logT_AGN": args.log_t_agn,
            "lmax": lmax,
            "lmax_untruncated": args.lmax_untruncated,
            "zmax": args.zmax,
            "version": args.version,
            "blind": args.blind,
            "reporting_ensemble": args.reporting_ensemble,
        }),
    )
    print(f"sigma_8 = {sigma8:.5f}, A_s = {a_s:.4e}")
    for i, t in enumerate(theta):
        print(f"{t:8.3f}  xip {xip[i]:12.5e}  xim {xim[i]:12.5e}  "
              f"band xi+ {xip[i] / xip_u[i]:6.4f} xi- {xim[i] / xim_u[i]:6.4f}  "
              f"flat/curved {xip_u[i] / xip_curved[i]:6.4f} "
              f"{xim_u[i] / xim_curved[i]:6.4f}")


if __name__ == "__main__":
    main()
