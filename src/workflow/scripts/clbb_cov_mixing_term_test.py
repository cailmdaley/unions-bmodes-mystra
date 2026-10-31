"""Isolate the E-to-B mixing term in the analytic C_l^BB bandpower covariance.

Replicates ``CosmologyValidation.calculate_pseudo_cl_eb_cov`` exactly -- same
mask (the mock catalogue's weighted galaxy-count map), nside, l range, powspace
binning, analytic noise bias, and the ``fiducial_input_inka = "coupled"``
branch -- but recomputes the covariance under three fiducial E-mode inputs
while sharing one NmtField / NmtWorkspace / NmtCovarianceWorkspace:

  a_theory  the chain's own input, [C_l^EE * pixwin^2, 0, 0, 0] + noise bias
  b_noEE    C_l^EE set to zero: noise bias only
  c_mockEE  the suite's own mean C_l^EE in place of theory

(a) is the gate: it must reproduce the stored covariance's sigma_BB. (b) says
how much of the analytic BB variance is E-to-B leakage rather than shape noise.
(c) separates a wrong input spectrum from a wrong coupling kernel.

The mock empirical sigma_BB (sample standard deviation over the suite's own BB
bandpowers) is carried alongside, so sigma_emp / sigma_ana is reported per band
for each variant.
"""

import argparse
import glob
import json

import healpy as hp
import numpy as np
import pymaster as nmt
from astropy.io import fits
from sp_validation.cosmology import get_cosmo, get_theo_c_ell


def build_bin(lmin, lmax, b_lmax, nbins, power):
    """The chain's powspace binning (sp_validation.pseudo_cl.make_namaster_bin)."""
    ells = np.arange(lmin, lmax + 1)
    edges = np.power(
        np.linspace(lmin**power, lmax**power, nbins + 1), 1.0 / power
    )
    bpws = np.digitize(ells.astype(float), edges) - 1
    bpws[0] = 0
    bpws[-1] = nbins - 1
    return nmt.NmtBin(ells=ells, bpws=bpws, lmax=b_lmax)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--mock-catalogs",
        required=True,
        help="glob over the realisation catalogues (the declared family input)",
    )
    ap.add_argument(
        "--realisation",
        default="00001",
        help="which realisation to build the mask and noise bias from",
    )
    ap.add_argument("--redshift-distribution", required=True)
    ap.add_argument("--cl-files", required=True, help="glob for the suite's C_l .npy")
    ap.add_argument("--cl-cov-mockmatched", required=True, help="stored cov, the gate")
    ap.add_argument("--nside", type=int, default=1024)
    ap.add_argument("--nbins", type=int, default=32)
    ap.add_argument("--power", type=float, default=0.5)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    nside = a.nside
    lmin, lmax, b_lmax = 8, 2 * nside, 2 * nside - 1
    b = build_bin(lmin, lmax, b_lmax, a.nbins, a.power)
    ell_eff = b.get_effective_ells()

    # ---- theory C_l^EE, on the chain's own (off-by-one) ell grid -------------
    ell = np.arange(1, lmax + 1)
    pw = hp.pixwin(nside, lmax=lmax)[1 : len(ell) + 1]
    z, dndz = np.loadtxt(a.redshift_distribution, unpack=True)
    cl_ee_theory = (
        get_theo_c_ell(ell=ell, z=z, nz=dndz, backend="ccl", cosmo=get_cosmo()) * pw**2
    )

    # ---- mask and noise bias from the mock catalogue ------------------------
    matches = [p for p in sorted(glob.glob(a.mock_catalogs)) if a.realisation in p]
    if len(matches) != 1:
        raise SystemExit(
            f"realisation {a.realisation} matched {len(matches)} of "
            f"{a.mock_catalogs}: {matches}"
        )
    mock_catalog = matches[0]
    print(f"mask and noise bias from {mock_catalog}", flush=True)

    with fits.open(mock_catalog, memmap=True) as hdul:
        cat = hdul["SOURCE_CATALOGUE"].data if "SOURCE_CATALOGUE" in hdul else hdul[1].data
        ra = np.array(cat["RA"], float)
        dec = np.array(cat["Dec"], float)
        e1 = np.array(cat["e1"], float)
        e2 = np.array(cat["e2"], float)
        w = np.array(cat["w"], float)

    pix = hp.ang2pix(nside, (90.0 - dec) * np.pi / 180.0, ra * np.pi / 180.0)
    unique_pix, _, idx_rep = np.unique(pix, return_index=True, return_inverse=True)
    n_gal = np.zeros(hp.nside2npix(nside))
    n_gal[unique_pix] = np.bincount(idx_rep, weights=w)

    variance_map = np.zeros(hp.nside2npix(nside))
    variance_map[unique_pix] = np.bincount(idx_rep, weights=(e1**2 + e2**2) / 2 * w**2)
    noise_bias = hp.nside2pixarea(nside) * np.mean(variance_map)

    zero = np.zeros(hp.nside2npix(nside))
    f = nmt.NmtField(mask=n_gal, maps=[zero, zero], lmax=b_lmax)
    wsp = nmt.NmtWorkspace.from_fields(f, f, b)
    cw = nmt.NmtCovarianceWorkspace.from_fields(f, f, f, f)

    # Noise bias: flat in EE and BB, decoupled, unbinned, filled below lmin --
    # exactly the chain's round trip.
    nb = np.zeros((4, lmax))
    nb[0, :] = noise_bias
    nb[3, :] = noise_bias
    nb = b.unbin_cell(wsp.decouple_cell(nb))
    lowest_ell = b.get_ell_list(0)[0]
    nb[:, :lowest_ell] = nb[:, [lowest_ell]]

    coupling = np.reshape(wsp.get_coupling_matrix(), (4, lmax, 4, lmax), order="F")
    mean_mask_sq = np.mean(n_gal**2)
    n_bands = b.get_n_bands()

    # ---- the suite's empirical spectra --------------------------------------
    files = sorted(glob.glob(a.cl_files))
    if not files:
        raise SystemExit(f"no C_l files matched {a.cl_files}")
    stack = np.array([np.load(p) for p in files])
    sigma_emp = stack[:, 4, :].std(axis=0, ddof=1)
    mock_ee_mean = stack[:, 1, :].mean(axis=0)
    if not np.allclose(stack[0, 0, :], ell_eff):
        raise SystemExit("suite C_l are not on the reporting band grid")

    # The mock mean EE lives on bands; spread it back over multipoles so it can
    # stand in for the theory spectrum on the same ell grid.
    mock_ee_ell = b.unbin_cell(np.array([mock_ee_mean, mock_ee_mean * 0,
                                         mock_ee_mean * 0, mock_ee_mean * 0]))[0]
    mock_ee_ell[:lowest_ell] = mock_ee_ell[lowest_ell]
    # chain's ell grid is arange(1, lmax+1); unbin_cell's is arange(0, lmax)
    mock_ee_on_chain_grid = np.concatenate([mock_ee_ell[1:], mock_ee_ell[-1:]])

    variants = {
        "a_theory": cl_ee_theory,
        "b_noEE": np.zeros_like(cl_ee_theory),
        "c_mockEE": mock_ee_on_chain_grid,
    }

    results = {}
    for name, cl_ee in variants.items():
        fid = np.array([cl_ee, 0.0 * cl_ee, 0.0 * cl_ee, 0.0 * cl_ee]) + nb
        fid = np.tensordot(coupling, fid) / mean_mask_sq
        cov = nmt.gaussian_covariance(
            cw, 2, 2, 2, 2, fid, fid, fid, fid, wsp, wb=wsp
        ).reshape([n_bands, 4, n_bands, 4])
        results[name] = np.sqrt(np.diag(cov[:, 3, :, 3]))
        print(f"{name}: sigma_BB[:4] = {results[name][:4]}", flush=True)

    with fits.open(a.cl_cov_mockmatched) as hdul:
        sigma_stored = np.sqrt(np.diag(hdul["COVAR_BB_BB"].data))

    gate = float(np.max(np.abs(results["a_theory"] / sigma_stored - 1.0)))
    print(f"gate: max |a_theory/stored - 1| = {gate:.4e}", flush=True)

    out = {
        "ell_eff": ell_eff.tolist(),
        "reporting_cut": [300.0, 1600.0],
        "n_realisations": len(files),
        "noise_bias": float(noise_bias),
        "gate_max_frac_deviation_from_stored": gate,
        "sigma_BB_stored": sigma_stored.tolist(),
        "sigma_BB_empirical": sigma_emp.tolist(),
        "variants": {
            name: {
                "sigma_BB": s.tolist(),
                "ratio_emp_over_ana": (sigma_emp / s).tolist(),
            }
            for name, s in results.items()
        },
        "sources": {
            "mock_catalog": mock_catalog,
            "realisation": a.realisation,
            "redshift_distribution": a.redshift_distribution,
            "cl_files": a.cl_files,
            "cl_cov_mockmatched": a.cl_cov_mockmatched,
        },
    }
    with open(a.out_file, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {a.out_file}")


if __name__ == "__main__":
    main()
