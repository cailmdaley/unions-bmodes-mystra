"""The values behind the ASTRA decision options.

A recipe renders `{decisions.<id>}` as an option id. Each table here maps a
decision's option ids to the values the scripts compute with, and is the
argparse `type=` for that decision's flag:

    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)

Excluded options stay in the tables, so a table always covers its decision.
tools/test_decision_values.py keeps the keys equal to the spec's option ids.

`SURVEY_PROPERTIES` are facts about the catalogs from Paper I, keyed by
`fiducial_catalog`, their areas read from the `survey_areas` input; they cover every catalog version CosmoCov runs on, and a version
and its leakage-corrected twin share one entry, read by `survey_properties`.
`NZ_FILE` is the fiducial catalog's n(z) per `fiducial_blind`, which every catalog
version uses.

`MOCK_SUITE` is how each GLASS suite samples its galaxies.

`REFERENCE_COVARIANCE` is the survey scale CosmoCov builds each reference
covariance at; the NaMaster reference covariance takes its window and noise from
a catalog summary instead.

Shell recipes read a value through the CLI, one key per decision it depends on:

    python decision_values.py reporting_binning bins_20
    python decision_values.py nz_file blind_a
"""

import sys
from pathlib import Path
from typing import NamedTuple

import yaml


class Grid(NamedTuple):
    """Logarithmic angular bins; separations in arcmin."""

    min_sep: float
    max_sep: float
    nbins: int


class MockSuite(NamedTuple):
    """How a GLASS suite samples its galaxies."""

    footprint: str | None  # catalog version whose galaxies set the mask; None: the fiducial catalog
    n_arcmin2: float | None  # galaxy density; None: matched to the fiducial catalog
    sigma_e: float | None  # per-component intrinsic dispersion; None: the fiducial catalog's


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

# The fiducial catalog's SOM n(z), one per blind. Every catalog version uses it:
# the n(z) enters theory and covariances only, and the B-mode covariance of a null
# signal does not depend on it.
NZ_PATH = "data/nz/nz_SP_v1.4.6.3_{blind}.txt"
NZ_FILE = {option: NZ_PATH.format(blind=blind) for option, blind in FIDUCIAL_BLIND.items()}

# Each catalog's survey area in deg², the `survey_areas` input.
SURVEY_AREA = yaml.safe_load((Path(__file__).resolve().parents[3] / "data/survey_areas.yaml").read_text())

SURVEY_PROPERTIES = {
    "size_cut": SurveyProperties(SURVEY_AREA["SP_v1.4.6.3"], 4.957279270321334, 0.3783979556382064),
    "initial": SurveyProperties(SURVEY_AREA["SP_v1.4.5"], 6.477062501504287, 0.39740419232655594),
    "masked": SurveyProperties(SURVEY_AREA["SP_v1.4.8"], 4.983323636055024, 0.3797397825805285),
    "relaxed_flags": SurveyProperties(SURVEY_AREA["SP_v1.4.11.3"], 6.2581259680050625, 0.380499698062307),
}


def survey_properties(version):
    """The survey scale of a catalog version, with or without `_leak_corr`."""
    catalogs = {cat_version: option for option, cat_version in FIDUCIAL_CATALOG.items()}
    return SURVEY_PROPERTIES[catalogs[version.removesuffix("_leak_corr")]]


# The survey scale behind each reference covariance. The mock-matched entry keeps
# the data's area and sets the suite's noise: the mocks' galaxy count over that
# area, and their total shape dispersion.
REFERENCE_COVARIANCE = {
    "mockmatched": SurveyProperties(SURVEY_PROPERTIES["size_cut"].area, 5.7930820, 0.3792561),
    "data_covariance": SURVEY_PROPERTIES["size_cut"],
}

# The published suite's sampling, and the data-matched suite's, which
# glass_mock_sampling.py derives from the fiducial catalog.
MOCK_SUITE = {
    "v1_4_6_3_v2": MockSuite("1.4.6", 6.0905, 0.2684),
    "data_matched": MockSuite(None, None, None),
    "v1_4_6": MockSuite(None, 6.0905, 0.2684),
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
    "survey_properties": SURVEY_PROPERTIES,
    "reference_covariance": REFERENCE_COVARIANCE,
    "mock_suite": MOCK_SUITE,
    "reporting_binning": REPORTING_BINNING,
    "integration_grid": INTEGRATION_GRID,
    "patch_count": PATCH_COUNT,
    "mode_count": MODE_COUNT,
    "harmonic_validation_binning": HARMONIC_VALIDATION_BINNING,
    "fiducial_scale_cuts": FIDUCIAL_ELL_CUT,
}

if __name__ == "__main__":
    value = TABLES[sys.argv[1]]
    for option in sys.argv[2:]:
        value = value[option]
    print(value)
