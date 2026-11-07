"""COSEBI covariance Tᵀ C_ξ T at the fiducial cuts and over the full range.

One npz, keyed like the modes it belongs to: ``{fiducial,full}_cov``, the
E_n ⊕ B_n matrix, with its ``_sigma_E``, ``_sigma_B`` and ``_scale_cut``. Only the
bin centres of the ξ_± file enter.
"""

import argparse

import numpy as np
import treecorr
from cosebis_transform import cosebis_covariance
from decision_values import INTEGRATION_GRID, MODE_COUNT, REPORTING_BINNING
from scale_cuts import add_scale_cut_convention_argument, resolve_fiducial_scale_cuts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--version", required=True)
    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)
    ap.add_argument("--integration-grid", type=INTEGRATION_GRID.__getitem__, required=True)
    ap.add_argument("--mode-count", type=MODE_COUNT.__getitem__, required=True)
    ap.add_argument("--xi-integration", required=True, help="1000-bin TreeCorr ξ_± dump")
    ap.add_argument("--cov-integration", required=True, help="1000-bin ξ_± covariance")
    ap.add_argument("--out-file", required=True, help="Output npz (lc {output})")
    add_scale_cut_convention_argument(ap)
    a = ap.parse_args(argv)

    reporting, nmodes = a.reporting_binning, a.mode_count
    scale_cuts = {
        "fiducial": resolve_fiducial_scale_cuts(reporting, a.scale_cut_convention),
        "full": (reporting.min_sep, reporting.max_sep),
    }
    min_sep, max_sep, nbins = a.integration_grid
    gg = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units="arcmin")
    gg.read(a.xi_integration)
    cov_xipm = np.loadtxt(a.cov_integration)

    payload = {"nmodes": nmodes, "version": a.version}
    for key, scale_cut in scale_cuts.items():
        cov = cosebis_covariance(gg, cov_xipm, nmodes, scale_cut)
        payload[f"{key}_cov"] = cov
        payload[f"{key}_sigma_E"] = np.sqrt(np.diag(cov[:nmodes, :nmodes]))
        payload[f"{key}_sigma_B"] = np.sqrt(np.diag(cov[nmodes:, nmodes:]))
        payload[f"{key}_scale_cut"] = np.array(scale_cut)
    with open(a.out_file, "wb") as f:
        np.savez(f, **payload)
    print(f"Saved COSEBI covariance to {a.out_file}")


if __name__ == "__main__":
    main()
