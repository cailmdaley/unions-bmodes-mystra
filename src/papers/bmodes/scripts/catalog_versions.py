"""The catalog versions a Paper II analysis compares, and their labels.

`fiducial_catalog` and `leakage_correction` select the primary shear catalogue.
The comparison roster is fixed in every universe; tools/test_decision_values.py
keeps it equal to the sweep matrices. A scan slot keeps its role when the
primary catalogue is raw: the fixed slots stay, and `scan_version` resolves
which catalogue each one holds.

As a CLI it prints the fiducial version, or the catalogue of one scan slot,
which is how recipes name a version:

    python catalog_versions.py size_cut applied
    python catalog_versions.py size_cut applied --scan-version SP_v1.4.8 [--fine-grid]
    python catalog_versions.py size_cut applied --nonfiducial
    python catalog_versions.py size_cut applied --comparison
"""

import argparse
from dataclasses import dataclass

from decision_values import FIDUCIAL_CATALOG

VERSIONS = (
    "SP_v1.4.5", "SP_v1.4.5_leak_corr",
    "SP_v1.4.6.3", "SP_v1.4.6.3_leak_corr",
    "SP_v1.4.6_ecut07", "SP_v1.4.6_ecut07_leak_corr",
    "SP_v1.4.8", "SP_v1.4.8_leak_corr",
    "SP_v1.4.11.3", "SP_v1.4.11.3_leak_corr",
    "SP_v1.4.11.3_ecut07", "SP_v1.4.11.3_ecut07_leak_corr",
)

VERSION_LABELS = {
    "SP_v1.4.5": "Initial",
    "SP_v1.4.6.3": "Size cuts",
    "SP_v1.4.8": "Masking",
    "SP_v1.4.11.3": "Relaxed flags",
}


@dataclass(frozen=True)
class Catalogs:
    """The catalog versions under one fiducial_catalog and leakage_correction."""

    catalog: str
    leakage: str

    @property
    def corrected(self):
        return self.leakage != "uncorrected"

    @property
    def suffix(self):
        return "_leak_corr" if self.corrected else ""

    @property
    def fiducial(self):
        return FIDUCIAL_CATALOG[self.catalog] + self.suffix

    @property
    def labels(self):
        return {version + self.suffix: label for version, label in VERSION_LABELS.items()}

    @property
    def comparison(self):
        """Non-ecut versions on the primary leakage branch."""
        return [v for v in VERSIONS if "_ecut" not in v and ("_leak_corr" in v) == self.corrected]

    @property
    def nonfiducial(self):
        return [v for v in VERSIONS if "_ecut" not in v and v != self.fiducial]

    @property
    def scan_slots(self):
        """The fixed non-fiducial scan slots, which are the corrected branch's."""
        return Catalogs(self.catalog, "applied").nonfiducial

    def scan_version(self, version, fine=False):
        """The catalogue a scan slot holds; the primary variant is produced separately."""
        if self.corrected:
            return version
        if fine:
            return version.removesuffix("_leak_corr")
        return self.fiducial + "_leak_corr" if version == self.fiducial else version


def add_catalog_arguments(parser):
    parser.add_argument("--fiducial-catalog", required=True, choices=FIDUCIAL_CATALOG)
    parser.add_argument("--leakage-correction", required=True, choices=["applied", "uncorrected"])


def catalogs_from(args):
    return Catalogs(args.fiducial_catalog, args.leakage_correction)


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(description="Print a Paper II catalog version.")
    ap.add_argument("catalog", choices=FIDUCIAL_CATALOG)
    ap.add_argument("leakage", choices=["applied", "uncorrected"])
    ap.add_argument("--scan-version", help="Resolve a fixed scan slot's catalogue")
    ap.add_argument("--fine-grid", action="store_true")
    ap.add_argument("--nonfiducial", action="store_true", help="List the non-fiducial versions")
    ap.add_argument("--comparison", action="store_true",
                    help="List the non-fiducial versions of the comparison figures")
    a = ap.parse_args(argv)
    catalogs = Catalogs(a.catalog, a.leakage)
    if a.nonfiducial:
        print("\n".join(catalogs.nonfiducial))
    elif a.comparison:
        print("\n".join(v for v in catalogs.comparison if v != catalogs.fiducial))
    elif a.scan_version:
        print(catalogs.scan_version(a.scan_version, a.fine_grid))
    else:
        print(catalogs.fiducial)


if __name__ == "__main__":
    _from_cli()
