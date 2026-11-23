"""The GLASS realizations: shell spectra, sampling and catalogs.

The configuration is the one the v1.4.6.3_v2 suite was published with
(sp_validation 4ec270d, ``glass_mock/make_unions_glass_sim.py``, run at nside
4096): Planck18
densities with n_s = 0.9665, m_ν = 0.06 eV taken from the CDM density, A_s
rescaled to σ8 = 0.8102, CAMB HMcode ``mead2020_feedback`` with log T_AGN = 7.8;
linear-window shells every 200 Mpc to z = 3; non-Limber matter Cℓ to ℓ = 4096
passed to the lognormal generator (ncorr = 3) without pixel-window
discretisation; one source population with the n(z) cut at z ≤ 3 and scaled to the
suite's density without renormalising; linear galaxy bias 1.2; Gaussian
photo-z with σ_0 = 0.03; Gaussian intrinsic ellipticities; unit weights.
Realization N is seeded with 1 000 000 + N.

The suite's density and per-component σ_e are its sampling
(``glass_mock_sampling.py``): the published 6.0905 arcmin⁻² and 0.2684 for
``v1_4_6_3_v2``; for ``data_matched``, the density at which a realization's
expected count on the footprint equals the fiducial catalog's weight-effective
count (Σw)²/Σw², and the catalog's σ_e.

The footprint is the nside-4096 occupancy of the suite's footprint catalog, read
from the catalog summary (``catalog_summary.py``).

    glass_mock_generate.py --number N --summary SUMMARY --suite SUITE --fiducial-catalog CAT \
        --nz NZ --shell-cls CLS --sampling S --out-file CAT

``glass_mock_shell_cls.py`` and ``glass_mock_sampling.py`` make the shell
spectra and the sampling from the functions here.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

NSIDE = 4096
LMAX = NSIDE
DX = 200.0
ZMAX = 3.0
BIAS = 1.2
PHZ_SIGMA_0 = 0.03
SEED_OFFSET = 1_000_000

H, OM, OB = 0.6766, 0.30966, 0.04897  # astropy Planck18
NS, SIGMA8, MNU, LOG_T_AGN, AS_INIT, KMAX = 0.9665, 0.8102, 0.06, 7.8, 2.1e-9, 20.0


def camb_params():
    """CAMB parameters with A_s rescaled so that σ8(z=0) = SIGMA8."""
    import camb

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=100 * H, omch2=(OM - OB) * H**2, ombh2=OB * H**2, mnu=MNU)
    oc = OM - OB - pars.omeganu
    pars.set_cosmology(H0=100 * H, omch2=oc * H**2, ombh2=OB * H**2, mnu=MNU)
    pars.InitPower.set_params(As=AS_INIT, ns=NS)
    pars.WantTransfer = True
    pars.NonLinear = camb.model.NonLinear_both
    pars.NonLinearModel.set_params(halofit_version="mead2020_feedback", HMCode_logT_AGN=LOG_T_AGN)
    pars.set_matter_power(nonlinear=True, kmax=KMAX)
    sigma8 = camb.get_results(pars).get_sigma8_0()
    power = camb.InitialPowerLaw()
    power.set_params(As=AS_INIT * (SIGMA8 / sigma8) ** 2, ns=NS)
    pars.InitPower = power
    return pars


def cosmology(pars):
    import camb
    from cosmology.compat.camb import Cosmology

    return Cosmology(camb.get_background(pars))


def shells_for(pars):
    import glass

    return glass.linear_windows(glass.distance_grid(cosmology(pars), 0.0, ZMAX, dx=DX))


def cut_nz(nz_file):
    """The n(z) the shells receive: cut at ZMAX, not renormalised."""
    nz = np.loadtxt(nz_file)
    nz = nz[nz[:, 0] <= ZMAX]
    return nz[:, 0], nz[:, 1]


def footprint_key(suite, fiducial_catalog):
    """The catalog-summary key of the catalog whose galaxies set a suite's
    footprint (and, when data-matched, its sampling): the suite's own, else the
    fiducial one."""
    from decision_values import FIDUCIAL_CATALOG, MOCK_SUITE

    footprint = MOCK_SUITE[suite].footprint
    return f"SP_v{footprint}" if footprint else FIDUCIAL_CATALOG[fiducial_catalog]


def footprint_mask(summary, suite, fiducial_catalog):
    """The nside-4096 pixels the suite samples in, as a float64 0/1 map."""
    from catalog_summary import footprint, load

    return footprint(load(summary, footprint_key(suite, fiducial_catalog)))


def sampling(suite, fiducial_catalog, summary, nz_file, out_file):
    """The suite's density and per-component σ_e, as JSON.

    A suite with published values keeps them. ``data_matched`` sets the density
    so that the expected count on the footprint, density × footprint area × the
    cut n(z) the shells receive, equals the fiducial catalog's (Σw)²/Σw², and
    σ_e² = Σw²(e1² + e2²) / (2 Σw²) on the raw ellipticities, the definition
    behind the data covariance's shape noise.
    """
    from decision_values import MOCK_SUITE

    published = MOCK_SUITE[suite]
    if published.n_arcmin2 is not None:
        out = {"n_arcmin2": published.n_arcmin2, "sigma_e": published.sigma_e}
    else:
        import glass
        from catalog_summary import load, survey

        entry = load(summary, footprint_key(suite, fiducial_catalog))
        area, n_eff, sigma_e = survey(entry)
        ngal = glass.partition(*cut_nz(nz_file), shells_for(camb_params()))
        out = {
            "n_arcmin2": float(n_eff / (area * 3600 * np.sum(ngal))),
            "sigma_e": sigma_e,
            "n_eff_catalog": n_eff,
            "galaxies_catalog": int(entry["n"]),
            "footprint_area_deg2": area,
        }
    with open(out_file, "w") as f:
        json.dump(out, f, indent=1)
    print(out)


def shell_cls(out_file):
    import glass.ext.camb

    pars = camb_params()
    cls = glass.ext.camb.matter_cls(pars, LMAX, shells_for(pars))
    np.save(out_file, np.asarray(cls))


def realization(number, mask, nz_file, cls_file, sampling_file, out_file):
    import fitsio
    import glass

    t0 = time.time()
    with open(sampling_file) as f:
        sample = json.load(f)
    n_arcmin2, sigma_e = sample["n_arcmin2"], sample["sigma_e"]
    pars = camb_params()
    shells = shells_for(pars)
    cls = list(np.load(cls_file))
    rng = np.random.default_rng(SEED_OFFSET + number)

    fields = glass.lognormal_fields(shells)
    gls = glass.solve_gaussian_spectra(fields, cls)
    matter = glass.generate(fields, gls, NSIDE, ncorr=3, rng=rng)
    convergence = glass.MultiPlaneConvergence(cosmology(pars))

    z, dndz = cut_nz(nz_file)
    dndz = dndz * n_arcmin2
    ngal = glass.partition(z, dndz, shells)
    zbins = glass.equal_dens_zbins(z, dndz, nbins=1)

    out_file = Path(out_file)
    part = out_file.with_name(out_file.name + ".part")
    fits = fitsio.FITS(str(part), "rw", clobber=True)
    fits.write(None)
    fits.create_table_hdu(
        names=["RA", "Dec", "e1", "e2", "w", "n1", "n2", "TOM_BIN_ID", "TRUE_Z", "PHOTO_Z"],
        formats=["D", "D", "E", "E", "D", "E", "E", "J", "D", "D"],
        extname="SOURCE_CATALOGUE",
    )
    hdu = fits["SOURCE_CATALOGUE"]
    hdu.write_key("GLASSVER", glass.__version__)
    hdu.write_key("SEED", SEED_OFFSET + number)
    hdu.write_key("MOCKNUM", number)
    hdu.write_key("NEFF", n_arcmin2)
    hdu.write_key("SIGMA_E", sigma_e)
    dtype = hdu.get_rec_dtype()[0]

    total = 0
    for i, delta in enumerate(matter):
        convergence.add_window(delta, shells[i])
        kappa = convergence.kappa
        g1, g2 = glass.shear_from_convergence(kappa)
        zeros = np.zeros_like(kappa)
        for lon, lat, count in glass.points.positions_from_delta(ngal[i], delta, BIAS, mask, rng=rng):
            total += count
            gal_z = glass.redshifts(count, shells[i], rng=rng)
            gal_phz = glass.gaussian_phz(gal_z, PHZ_SIGMA_0, rng=rng)
            ellip = glass.ellipticity_intnorm(count, sigma_e, rng=rng, xp=np)
            she = glass.galaxy_shear(lon, lat, ellip, kappa, g1, g2)
            noise = glass.galaxy_shear(lon, lat, ellip, zeros, zeros, zeros)
            cat = np.empty(count, dtype=dtype)
            cat["RA"], cat["Dec"] = lon, lat
            cat["e1"], cat["e2"] = she.real, -she.imag
            cat["w"] = 1.0
            cat["n1"], cat["n2"] = noise.real, noise.imag
            cat["TOM_BIN_ID"] = np.digitize(gal_phz, np.unique(zbins)) - 1
            cat["TRUE_Z"], cat["PHOTO_Z"] = gal_z, gal_phz
            hdu.append(cat)
    fits.close()
    part.rename(out_file)
    print(f"{out_file}: {total} galaxies in {time.time() - t0:.0f}s", flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--number", type=int, required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--suite", required=True)
    p.add_argument("--fiducial-catalog", required=True)
    p.add_argument("--nz", required=True)
    p.add_argument("--shell-cls", required=True)
    p.add_argument("--sampling", required=True)
    p.add_argument("--out-file", required=True)
    a = p.parse_args(argv)
    mask = footprint_mask(a.summary, a.suite, a.fiducial_catalog)
    realization(a.number, mask, a.nz, a.shell_cls, a.sampling, a.out_file)


if __name__ == "__main__":
    main()
