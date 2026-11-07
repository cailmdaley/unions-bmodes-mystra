"""Pin the option-to-value tables to the spec and to the published values.

`definition_version` does not hash `decision_values.py`, so editing a table
stales nothing; these tests are what hold it. Run on the laptop with PyYAML:

    python tools/test_decision_values.py
"""
import ast
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src/papers/bmodes/scripts"), str(ROOT / "src/workflow/scripts")]

import catalog_versions
import decision_values as dv
import generate_cosmocov_ini
import mock_vs_data_noise_table
import scale_cuts

ROOT_DECISIONS = yaml.safe_load((ROOT / "astra.yaml").read_text())["decisions"]
CAT_CONFIG = yaml.safe_load((ROOT / "src/cosmo_val/cat_config.yaml").read_text())
REAL_SPACE_DECISIONS = yaml.safe_load((ROOT / "analyses/data_vectors/astra.yaml").read_text())["decisions"]

# Tables keyed by the options of a decision other than their namesake.
KEYED_BY = {"nz_file": "fiducial_catalog", "survey_properties": "fiducial_catalog",
            "reference_catalog": "reference_covariance"}

FIDUCIAL_NZ = "/n17data/sguerrini/UNIONS/WL/nz/v1.4.6.3/nz_SP_v1.4.6.3_{}.txt"

# The values Paper II was computed with.
PUBLISHED = {
    "fiducial_blind": ("blind_a", "A"),
    "fiducial_catalog": ("size_cut", "SP_v1.4.6.3"),
    "nz_file": ("size_cut", {"blind_a": FIDUCIAL_NZ.format("A"), "blind_b": FIDUCIAL_NZ.format("B"),
                             "blind_c": FIDUCIAL_NZ.format("C")}),
    "survey_properties": ("size_cut", (2894.0303815287743, 4.957279270321334, 0.3783979556382064)),
    "reference_covariance": ("mockmatched", (2894.0303815287743, 5.793082, 0.3792561)),
    "reference_catalog": ("mockmatched", "SP_v1.4.6.3_glass_mock_v2"),
    "bb_blind_check_version": ("size_cut_v1463", "SP_v1.4.6.3_leak_corr"),
    "reporting_binning": ("bins_20", (1.0, 250.0, 20)),
    "integration_grid": ("bins_1000", (0.5, 300.0, 1000)),
    "patch_count": ("npatch_1", 1),
    "mode_count": ("n_max_20", 20),
    "harmonic_validation_binning": ("bins_96", 96),
    "mc_samples": ("n_2000", 2000),
    "fiducial_scale_cuts": ("adopted", (300, 1600)),
}


class Tables(unittest.TestCase):
    def test_keys_are_the_spec_option_ids(self):
        for name, table in dv.TABLES.items():
            decisions = REAL_SPACE_DECISIONS if name == "patch_count" else ROOT_DECISIONS
            self.assertEqual(set(table), set(decisions[KEYED_BY.get(name, name)]["options"]), name)

    def test_nz_files_cover_every_blind(self):
        for blinds in dv.NZ_FILE.values():
            self.assertEqual(set(blinds), set(ROOT_DECISIONS["fiducial_blind"]["options"]))

    def test_catalog_facts_match_cat_config(self):
        for option, version in dv.FIDUCIAL_CATALOG.items():
            entry = CAT_CONFIG[version]
            self.assertEqual(dv.NZ_FILE[option]["blind_a"], entry["shear"]["redshift_path"], option)
            cov_th = entry["cov_th"]
            self.assertEqual(dv.SURVEY_PROPERTIES[option], (cov_th["A"], cov_th["n_e"], cov_th["sigma_e"]), option)

    def test_every_cosmocov_version_has_survey_properties(self):
        for fiducial in dv.FIDUCIAL_CATALOG:
            for leakage in ("applied", "uncorrected"):
                catalogs = catalog_versions.Catalogs(fiducial, leakage)
                for version in [catalogs.fiducial, *map(catalogs.scan_version, catalogs.scan_slots)]:
                    cov_th = CAT_CONFIG[version.removesuffix("_leak_corr")]["cov_th"]
                    self.assertEqual(dv.survey_properties(version),
                                     (cov_th["A"], cov_th["n_e"], cov_th["sigma_e"]), version)

    def test_reference_catalogs_are_in_cat_config(self):
        for option, version in dv.REFERENCE_CATALOG.items():
            self.assertIn(version.removesuffix("_leak_corr"), CAT_CONFIG, option)

    def test_published_values(self):
        self.assertEqual(set(PUBLISHED), set(dv.TABLES))
        for name, (option, value) in PUBLISHED.items():
            self.assertEqual(dv.TABLES[name][option], value, name)

    def test_grids_carry_float_extents(self):
        # Tar member names embed these: `minsep=0.5_maxsep=300.0_nbins=1000`.
        for table in (dv.REPORTING_BINNING, dv.INTEGRATION_GRID):
            for grid in table.values():
                self.assertIsInstance(grid.min_sep, float)
                self.assertIsInstance(grid.max_sep, float)
                self.assertIsInstance(grid.nbins, int)

    def test_labels_quote_the_values(self):
        for name, table in (("reporting_binning", dv.REPORTING_BINNING),
                            ("integration_grid", dv.INTEGRATION_GRID)):
            for option, grid in table.items():
                label = ROOT_DECISIONS[name]["options"][option]["label"].replace(",", "")
                self.assertIn(str(grid.nbins), label, option)
        for name, table in (("mc_samples", dv.MC_SAMPLES), ("mode_count", dv.MODE_COUNT),
                            ("harmonic_validation_binning", dv.HARMONIC_VALIDATION_BINNING)):
            for option, value in table.items():
                self.assertIn(str(value), ROOT_DECISIONS[name]["options"][option]["label"], option)

    def test_scale_cut_edges_and_window(self):
        grid = dv.REPORTING_BINNING["bins_20"]
        edges = np.geomspace(1.0, 250.0, 21)
        self.assertEqual(scale_cuts.resolve_fiducial_scale_cuts(grid, "grid_bin_edges"),
                         (float(edges[9]), float(edges[16])))
        self.assertEqual(scale_cuts.resolve_fiducial_scale_cuts(grid, "rounded_arcmin"), (12.0, 83.0))
        for convention in scale_cuts.CONVENTIONS:
            self.assertEqual(scale_cuts.resolve_fiducial_bin_window(grid, convention), (9, 15))


class ReferenceNoise(unittest.TestCase):
    RECORDED_INI = ROOT / "analyses/covariances/results/baseline/cov_integration_g/covariance.ini"

    def ini(self, tmp, *flags):
        generate_cosmocov_ini.main([
            "--version", "SP_v1.4.6.3_leak_corr", "--blind", "blind_a",
            "--planck18-json", str(ROOT / "analyses/covariances/results/baseline/cosmology_snapshot/planck18.json"),
            "--min-sep", "0.5", "--max-sep", "300.0", "--nbins", "1000", "--gaussian", "g",
            "--mask-cls", "/n17data/cdaley/unions/code/sp_validation/cosmo_inference/data/mask/mask_cls_footprint_nside_4096_norm.txt",
            "--nz", dv.NZ_FILE["size_cut"]["blind_a"], "--out-ini", f"{tmp}/covariance.ini", *flags,
        ])
        return Path(tmp, "covariance.ini").read_text().splitlines()

    def test_data_option_is_the_recorded_ini(self):
        recorded = self.RECORDED_INI.read_text().splitlines()
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.ini(tmp), recorded)
            self.assertEqual(self.ini(tmp, "--reference-covariance", "data_covariance"), recorded)

    def test_mock_option_changes_only_the_noise(self):
        recorded = self.RECORDED_INI.read_text().splitlines()
        with tempfile.TemporaryDirectory() as tmp:
            changed = [new for new, old in zip(self.ini(tmp, "--reference-covariance", "mockmatched"), recorded)
                       if new != old]
        self.assertEqual(changed, ["sigma_e : 0.3792561", "source_n_gal : 5.793082", "lens_n_gal : 5.793082"])

    def test_noise_table_ratios(self):
        audit = dict(data_galaxies=61378891, data_area_deg2=2750.65, mock_galaxies=60423575, mock_area_deg2=2709.91)
        mock = mock_vs_data_noise_table.build_rows(
            dv.REFERENCE_COVARIANCE["data_covariance"], dv.REFERENCE_COVARIANCE["mockmatched"], audit)[1]
        self.assertEqual(f"{mock['noise_power_ratio_to_data']:.6g}", "0.85961")
        self.assertEqual(f"{mock['bb_variance_ratio_to_data']:.6g}", "0.738929")
        self.assertEqual(mock["galaxies_from_n_gal_convention"], 60355279)


class Roster(unittest.TestCase):
    def test_scan_matrices_stay_inside_the_roster(self):
        sidecars = list((ROOT / "analyses").rglob("*.sweep.yaml"))
        found = 0
        for path in sidecars:
            matrix = yaml.safe_load(path.read_text()).get("matrix") or {}
            for version in matrix.get("version") or []:
                if str(version).startswith("SP_v"):
                    found += 1
                    self.assertIn(version, catalog_versions.VERSIONS, path.name)
        self.assertGreater(found, 0)

    def test_xi_scan_matrix_is_the_scan_slots(self):
        matrix = yaml.safe_load(
            (ROOT / "analyses/data_vectors/analyses/xi_scan/xi_scan.sweep.yaml").read_text()
        )["matrix"]
        slots = catalog_versions.Catalogs("size_cut", "applied").scan_slots
        self.assertEqual(set(matrix["version"]), set(slots))


class Plumbing(unittest.TestCase):
    def test_config_cosebis_always_receives_the_integration_grid(self):
        source = ROOT / "src/papers/bmodes/scripts/harmonic_config_cosebis_comparison.py"
        calls = [node for node in ast.walk(ast.parse(source.read_text()))
                 if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "_compute_config_cosebis"]
        self.assertTrue(calls)
        for call in calls:
            self.assertEqual(call.args[-1].id, "integration", call.lineno)


class NoAnalysisConfig(unittest.TestCase):
    """The analysis has no config file: values arrive as flags or module constants."""

    SOURCES = sorted(
        p for pattern in ("src/papers/bmodes/scripts/*", "src/workflow/scripts/*")
        for p in ROOT.glob(pattern) if p.suffix in (".py", ".sh") and p.is_file() and not p.is_symlink()
    )
    SPECS = [ROOT / "astra.yaml", *(ROOT / "analyses").rglob("astra.yaml"),
             *(ROOT / "analyses").rglob("*.sweep.yaml")]

    def test_removed_names_stay_removed(self):
        banned = re.compile(r"selected_config|paper_config|sweep_versions|\bload_config|select_leakage"
                            r"|bmodes/config/config\.yaml|config/catalogs/")
        for path in [*self.SOURCES, *self.SPECS]:
            self.assertIsNone(banned.search(path.read_text()), path.name)

    def test_no_config_style_lookups(self):
        keys = {"plotting", "x_offsets", "version_alpha", "version_box", "cosebis_nbins", "version_labels"}
        for path in self.SOURCES:
            if path.suffix != ".py":
                continue
            for node in ast.walk(ast.parse(path.read_text())):
                if (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
                        and node.value.id in ("config", "plotting_config")):
                    self.fail(f"{path.name}:{node.lineno} indexes {node.value.id}")
                if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                        and node.func.attr == "get" and node.args
                        and isinstance(node.args[0], ast.Constant)
                        and node.args[0].value in keys):
                    self.fail(f"{path.name}:{node.lineno} reads a config key with .get")


if __name__ == "__main__":
    unittest.main()
