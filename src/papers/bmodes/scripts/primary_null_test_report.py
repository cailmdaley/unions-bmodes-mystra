"""Select the primary reported range and recalibrate the selected test family.

The full six-test, two-range evidence stays available as a reference product.
Dropping a COSEBI subset recalibrates the minimum from stored mock PTE columns;
it never remeasures or retransforms the mock catalogues.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "workflow/scripts"))
from mock_min_pte_calibration import _calibration_payload


SUBSETS = {
    "n_le_6_and_n_le_20": {"B_n<=6", "B_n<=20"},
    "n_le_6_only": {"B_n<=6"},
    "n_le_20_only": {"B_n<=20"},
}


def select_report(reference, mock_names, mock_ptes, range_option, subset_option):
    cut = "fiducial" if range_option == "adopted" else "full_range"
    ref = reference if cut == "fiducial" else reference["full_range"]
    keep = [i for i, name in enumerate(mock_names)
            if not name.startswith("B_n<=") or name in SUBSETS[subset_option]]
    names = [mock_names[i] for i in keep]
    observed = {name: ref["observed_ptes"][name] for name in names}
    payload = _calibration_payload(
        mock_ptes[:, keep], names, observed, tuple(ref["window_bins"]),
        tuple(ref["scale_cut_arcmin"]), ref["scale_cut_convention"], ref["global_alpha"],
    )
    payload["definition"] = (
        "Empirical fraction of clean mocks whose minimum PTE across the selected "
        "correlated tests is <= the observed minimum."
    )
    return dict(primary_range=cut, mode_subsets=subset_option, calibration=payload)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calibration", required=True, type=Path)
    ap.add_argument("--per-mock-statistics", required=True, type=Path)
    ap.add_argument("--fiducial-scale-cuts", required=True, choices=["adopted", "full_range"])
    ap.add_argument("--mode-subsets", required=True, choices=SUBSETS)
    ap.add_argument("--out-file", required=True, type=Path)
    a = ap.parse_args()
    reference = json.loads(a.calibration.read_text())
    with np.load(a.per_mock_statistics) as data:
        names = [str(name) for name in data["stat_names"]]
        ptes = data["pte" if a.fiducial_scale_cuts == "adopted" else "pte_full"]
    payload = select_report(reference, names, ptes, a.fiducial_scale_cuts, a.mode_subsets)
    a.out_file.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
