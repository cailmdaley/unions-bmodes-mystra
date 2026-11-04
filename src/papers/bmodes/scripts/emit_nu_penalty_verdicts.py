"""Emit the ν → ν−2 scale-cut penalty verdict table for the paper's PTE table.

Referee comment 1 asks what happens to the B-mode null-test verdicts if the two
fiducial scale-cut boundaries are treated as *fitted* parameters, i.e. if the
chi-square degrees of freedom are reduced by two. This generalizes the
single-number footnote reproduced by ``emit_pte_combined_dof_penalized.py`` to
every cell of the paper's PTE table: four leak-corrected catalog versions x six
statistics x two cuts (fiducial and full range).

Every PTE in the table is an analytic ``scipy.stats.chi2.sf(chi2, nu)``:

* xi_+^B / xi_-^B   nu = number of theta bins in the cut window
                    (``sp_validation.b_modes.calculate_eb_statistics``);
* xi_tot^B          nu = 2 x that, the joint xi_+/xi_- test;
* COSEBIs           nu = n_modes (6 or 20), independent of the angular cut;
* C_l^BB            nu = number of bandpowers in the ell window
                    (``plotting_utils.compute_chi2_pte``, dof = len(data)).

nu is therefore fully determined by the cut, so chi2_obs is recovered exactly by
inverting the stored PTE and the penalized PTE follows as sf(chi2_obs, nu-2). No
statistic is recomputed and no MC is re-run.

The fiducial pure-mode window is resolved by ``scale_cuts.resolve_fiducial_bin_window``
and the fiducial ell window by ``plotting_utils.ell_bin_index`` -- the same
resolvers the PTE evidence itself used, never bin means.

Emits a CSV with one row per (version, statistic, cut).
"""

import argparse
import csv
import json
import sys

import numpy as np
import yaml
from scipy import stats

from plotting_utils import ell_bin_index
from scale_cuts import add_scale_cut_convention_argument, resolve_fiducial_bin_window
from paper_config import load_config

N_FITTED_BOUNDARIES = 2
PASS_THRESHOLD = 0.05

# Table row order, matching the paper's table: display label -> evidence key.
VERSIONS = [
    ("Initial", "SP_v1.4.5_leak_corr"),
    ("Size_cuts", "SP_v1.4.6.3_leak_corr"),
    ("Masked", "SP_v1.4.8_leak_corr"),
    ("Relaxed_flags", "SP_v1.4.11.3_leak_corr"),
]

FULL_RANGE_THETA_BINS = 20  # pte_matrix[0, len(theta)-1] over the 20-bin grid


def config_space_rows(label, version, evidence, n_fid_bins):
    """(statistic, cut, nu, pte) rows read from the config-space PTE evidence."""
    spec = {
        "cosebis_n6": ("cosebis_stats", 6, 6),
        "cosebis_n20": ("cosebis_20_stats", 20, 20),
        "xip_B": ("xip_stats", n_fid_bins, FULL_RANGE_THETA_BINS),
        "xim_B": ("xim_stats", n_fid_bins, FULL_RANGE_THETA_BINS),
        "xitot_B": ("combined_stats", 2 * n_fid_bins, 2 * FULL_RANGE_THETA_BINS),
    }
    for statistic, (key, nu_fid, nu_full) in spec.items():
        for cut, nu in (("fiducial", nu_fid), ("full_range", nu_full)):
            yield statistic, cut, nu, float(evidence[key][f"pte_at_{cut}"])


def harmonic_rows(evidence):
    """(statistic, cut, nu, pte) rows for C_l^BB."""
    ell = np.asarray(evidence["ell"], dtype=float)
    ell_min, ell_max = evidence["fiducial_ell_range"]
    i_min, i_max = ell_bin_index(ell, ell_min, ell_max)
    yield "cl_BB", "fiducial", i_max - i_min + 1, float(evidence["pte_at_fiducial"])
    yield "cl_BB", "full_range", len(ell), float(evidence["pte_at_full_range"])


def penalized_rows(config_space, harmonic, n_fid_bins, uncorrected=False):
    """One row per (version × statistic × cut) with its ν−2 PTE.

    *config_space* / *harmonic* are the ``evidence.versions`` mappings of the
    configuration- and harmonic-space PTE evidence.
    """
    rows = []
    for label, version in VERSIONS:
        if uncorrected:
            version = version.removesuffix("_leak_corr")
        entries = list(
            config_space_rows(label, version, config_space[version], n_fid_bins)
        )
        entries += list(harmonic_rows(harmonic[version]))
        for statistic, cut, nu, pte in entries:
            chi2_obs = float(stats.chi2.isf(pte, nu))
            nu_pen = nu - N_FITTED_BOUNDARIES
            pte_pen = float(stats.chi2.sf(chi2_obs, nu_pen))
            passes = pte >= PASS_THRESHOLD
            passes_pen = pte_pen >= PASS_THRESHOLD
            rows.append(
                {
                    "version": label,
                    "version_full": version,
                    "cut": cut,
                    "statistic": statistic,
                    "nu": nu,
                    "nu_penalized": nu_pen,
                    "chi2": f"{chi2_obs:.6g}",
                    "pte": f"{pte:.6g}",
                    "pte_penalized": f"{pte_pen:.6g}",
                    "pte_decrease": f"{pte - pte_pen:.6g}",
                    "passes": int(passes),
                    "passes_penalized": int(passes_pen),
                    "verdict_change": int(passes != passes_pen),
                }
            )
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config-space-evidence", required=True)
    p.add_argument("--harmonic-space-evidence", required=True)
    p.add_argument("--config", required=True, help="bmodes config.yaml")
    add_scale_cut_convention_argument(p)
    p.add_argument("--out-file", required=True)
    args = p.parse_args()

    config = load_config(args.config)
    first_bin, last_bin = resolve_fiducial_bin_window(config, args.scale_cut_convention)
    n_fid_bins = last_bin + 1 - first_bin

    with open(args.config_space_evidence) as f:
        config_space = json.load(f)["evidence"]["versions"]
    with open(args.harmonic_space_evidence) as f:
        harmonic = json.load(f)["evidence"]["versions"]

    rows = penalized_rows(config_space, harmonic, n_fid_bins,
                          uncorrected=config.get("leakage_correction") == "uncorrected")

    with open(args.out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fiducial = [r for r in rows if r["cut"] == "fiducial"]
    changed = [r for r in rows if r["verdict_change"]]
    print(
        f"n_fid_bins={n_fid_bins} rows={len(rows)} "
        f"verdict_changes={len(changed)} "
        f"max_fiducial_pte_decrease={max(float(r['pte_decrease']) for r in fiducial):.4f}"
    )
    for r in changed:
        print(f"  CHANGED: {r['version']} {r['statistic']} {r['cut']} "
              f"{r['pte']} -> {r['pte_penalized']}", file=sys.stderr)


if __name__ == "__main__":
    main()
