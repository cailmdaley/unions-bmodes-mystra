"""Emit the degrees-of-freedom-penalized combined B-mode PTE as an ASTRA
`type: metric` output.

The paper footnote notes that if the two fitted scale-cut boundaries are
treated as fitted parameters (reducing the chi-square degrees of freedom by 2),
the minimum combined PTE across all statistics rises from the un-penalized
value to a penalized one, still above the 0.05 threshold.

The combined PTE is `scipy.stats.chi2.sf(chi2_obs, nu)` with nu = 2 * (number
of theta bins in the cut) for the joint xi_+ / xi_- test. The per-cut PTE grid
npz stores the PTE matrix but not the raw chi2/nu; however nu is fully
determined by the fiducial scale cut (12-83 arcmin -> bins [9, 15] -> 7 bins
-> nu = 14), so chi2_obs is recovered exactly by inverting the stored
un-penalized PTE, and the penalized PTE follows with nu-2. No MC re-run needed.

The window is resolved by `scale_cuts.resolve_fiducial_bin_window`, the single
definition shared with the PTE matrices and the mock calibration.
"""

import argparse
import json
import os

import numpy as np
import yaml
from scipy import stats

from scale_cuts import (
    add_scale_cut_convention_argument,
    resolve_fiducial_bin_window,
    resolve_fiducial_scale_cuts,
)
from paper_config import load_config

# The two scale-cut boundaries are the "fitted" parameters.
N_FITTED_BOUNDARIES = 2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--pte-per-cut-npz",
        required=True,
        help="<version>_A_pure_eb_ptes.npz with theta + pte_combined grid",
    )
    p.add_argument("--config", required=True, help="bmodes config.yaml")
    add_scale_cut_convention_argument(p)
    p.add_argument("--out-file", required=True, help="output metric JSON path")
    args = p.parse_args()

    config = load_config(args.config)
    cut_arcmin = resolve_fiducial_scale_cuts(config, args.scale_cut_convention)
    s, e = resolve_fiducial_bin_window(config, args.scale_cut_convention)

    d = np.load(args.pte_per_cut_npz)
    nbins = (e + 1) - s
    nu = 2 * nbins  # xi_+ and xi_- jointly

    pte_unpen = float(d["pte_combined"][s, e])
    chi2_obs = float(stats.chi2.isf(pte_unpen, nu))
    nu_pen = nu - N_FITTED_BOUNDARIES
    pte_pen = float(stats.chi2.sf(chi2_obs, nu_pen))

    # Metric payload with provenance fields alongside (extra keys are
    # ignored by metric readers).
    metric = {
        "value": round(pte_pen, 4),
        "uncertainty": None,
        "unit": None,
        "statistic": "combined_xitot_B",
        "cut": "fiducial",
        "scale_cut_arcmin": list(cut_arcmin),
        "scale_cut_convention": args.scale_cut_convention,
        "bin_indices": [s, e],
        "nu": nu,
        "nu_penalized": nu_pen,
        "n_fitted_boundaries": N_FITTED_BOUNDARIES,
        "pte_unpenalized": pte_unpen,
        "chi2_obs": chi2_obs,
        "pte_penalized": pte_pen,
    }
    with open(args.out_file, "w") as f:
        json.dump(metric, f, indent=2)

    print(
        f"nu={nu} chi2={chi2_obs:.4f} pte_unpen={pte_unpen:.5f} pte_penalized={pte_pen:.5f}"
    )


if __name__ == "__main__":
    main()
