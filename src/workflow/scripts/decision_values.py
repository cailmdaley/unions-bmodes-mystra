"""The values behind the ASTRA decision options.

A recipe renders `{decisions.<id>}` as an option id. Each table here maps a
decision's option ids to the values the scripts compute with, and is the
argparse `type=` for that decision's flag:

    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)

Excluded options stay in the tables, so a table always covers its decision.
tools/test_decision_values.py keeps the keys equal to the spec's option ids.

`NZ_FILE`, `MASK_CLS` and `SURVEY_PROPERTIES` are facts about the catalogs from
Paper I, keyed by `fiducial_catalog` and, for the n(z), by `fiducial_blind`.
`SURVEY_PROPERTIES` covers every catalog version CosmoCov runs on; a version and
its leakage-corrected twin share one entry, read by `survey_properties`.

Shell recipes read a value through the CLI, one key per decision it depends on:

    python decision_values.py mc_samples n_2000
    python decision_values.py nz_file size_cut blind_a
"""

import sys
from typing import NamedTuple


class Grid(NamedTuple):
    """Logarithmic angular bins; separations in arcmin."""

    min_sep: float
    max_sep: float
    nbins: int


class SurveyProperties(NamedTuple):
    """Survey scale of a catalog."""

    area: float  # deg²
    n_eff: float  # arcmin⁻²
    sigma_e_combined: float  # two-component shape dispersion


FIDUCIAL_BLIND = {"blind_a": "A", "blind_b": "B", "blind_c": "C"}

# The cat_config version of each catalog; `_leak_corr` is added per leakage_correction.
FIDUCIAL_CATALOG = {
    "size_cut": "SP_v1.4.6.3",
    "initial": "SP_v1.4.5",
    "masked": "SP_v1.4.8",
    "relaxed_flags": "SP_v1.4.11.3",
}

# Each catalog's n(z) file per blind. The masked and relaxed-flags catalogs use
# the pre-SOM n(z).
NZ_FILE = {
    catalog: {option: path.format(blind=blind) for option, blind in FIDUCIAL_BLIND.items()}
    for catalog, path in {
        "size_cut": "/n17data/sguerrini/UNIONS/WL/nz/v1.4.6.3/nz_SP_v1.4.6.3_{blind}.txt",
        "initial": "/n17data/sguerrini/UNIONS/WL/nz/v1.4.5/nz_SP_v1.4.5_{blind}.txt",
        "masked": "/n17data/mkilbing/astro/data/CFIS/v1.0/nz/dndz_SP_{blind}.txt",
        "relaxed_flags": "/n17data/mkilbing/astro/data/CFIS/v1.0/nz/dndz_SP_{blind}.txt",
    }.items()
}

# Footprint mask power spectrum (nside 4096, normalised); only the fiducial
# catalog has one.
MASK_CLS = {
    "size_cut": "/n17data/cdaley/unions/code/sp_validation/cosmo_inference/data/mask/mask_cls_footprint_nside_4096_norm.txt",
    "initial": None,
    "masked": None,
    "relaxed_flags": None,
}

SURVEY_PROPERTIES = {
    "size_cut": SurveyProperties(2894.0303815287743, 4.957279270321334, 0.3783979556382064),
    "initial": SurveyProperties(2894.0303815287743, 6.477062501504287, 0.39740419232655594),
    "masked": SurveyProperties(2517.021356326454, 4.983323636055024, 0.3797397825805285),
    "relaxed_flags": SurveyProperties(2894.0303815287743, 6.2581259680050625, 0.380499698062307),
}


def survey_properties(version):
    """The survey scale of a catalog version, with or without `_leak_corr`."""
    catalogs = {cat_version: option for option, cat_version in FIDUCIAL_CATALOG.items()}
    return SURVEY_PROPERTIES[catalogs[version.removesuffix("_leak_corr")]]


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
    "nz_file": NZ_FILE,
    "mask_cls": MASK_CLS,
    "survey_properties": SURVEY_PROPERTIES,
    "reporting_binning": REPORTING_BINNING,
    "integration_grid": INTEGRATION_GRID,
    "patch_count": PATCH_COUNT,
    "mode_count": MODE_COUNT,
    "harmonic_validation_binning": HARMONIC_VALIDATION_BINNING,
    "mc_samples": MC_SAMPLES,
    "fiducial_scale_cuts": FIDUCIAL_ELL_CUT,
}

if __name__ == "__main__":
    value = TABLES[sys.argv[1]]
    for option in sys.argv[2:]:
        value = value[option]
    print(value)
