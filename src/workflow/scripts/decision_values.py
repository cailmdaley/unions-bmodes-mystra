"""The values behind the ASTRA decision options.

A recipe renders `{decisions.<id>}` as an option id. Each table here maps a
decision's option ids to the values the scripts compute with, and is the
argparse `type=` for that decision's flag:

    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)

Excluded options stay in the tables, so a table always covers its decision.
tools/test_decision_values.py keeps the keys equal to the spec's option ids.

Shell recipes read a scalar through the CLI:

    python decision_values.py mc_samples n_2000
"""

import sys
from typing import NamedTuple


class Grid(NamedTuple):
    """Logarithmic angular bins; separations in arcmin."""

    min_sep: float
    max_sep: float
    nbins: int


FIDUCIAL_BLIND = {"blind_a": "A", "blind_b": "B", "blind_c": "C"}

# The cat_config version of each catalog; `_leak_corr` is added per leakage_correction.
FIDUCIAL_CATALOG = {
    "size_cut": "SP_v1.4.6.3",
    "initial": "SP_v1.4.5",
    "masked": "SP_v1.4.8",
    "relaxed_flags": "SP_v1.4.11.3",
}

REPORTING_BINNING = {
    "bins_20": Grid(1.0, 250.0, 20),
    "bins_40": Grid(1.0, 250.0, 40),
}

INTEGRATION_GRID = {
    "bins_1000": Grid(0.5, 300.0, 1000),
    "bins_10000": Grid(0.5, 300.0, 10000),
}

PATCH_COUNT = {"npatch_1": 1, "npatch_120": 120}

MODE_COUNT = {"n_max_20": 20, "n_max_15": 15, "n_max_30": 30}

HARMONIC_VALIDATION_BINNING = {"bins_96": 96, "bins_128": 128, "bins_32": 32}

MC_SAMPLES = {"n_2000": 2000, "n_1000": 1000, "n_5000": 5000}

# Harmonic ℓ window of each fiducial_scale_cuts option. The B-mode subtraction
# option is a treatment, not a window.
FIDUCIAL_ELL_CUT = {
    "adopted": (300, 1600),
    "full_range": (8, 2048),
    "bmode_bias_subtraction": None,
}

TABLES = {
    "fiducial_blind": FIDUCIAL_BLIND,
    "fiducial_catalog": FIDUCIAL_CATALOG,
    "reporting_binning": REPORTING_BINNING,
    "integration_grid": INTEGRATION_GRID,
    "patch_count": PATCH_COUNT,
    "mode_count": MODE_COUNT,
    "harmonic_validation_binning": HARMONIC_VALIDATION_BINNING,
    "mc_samples": MC_SAMPLES,
    "fiducial_scale_cuts": FIDUCIAL_ELL_CUT,
}

if __name__ == "__main__":
    print(TABLES[sys.argv[1]][sys.argv[2]])
