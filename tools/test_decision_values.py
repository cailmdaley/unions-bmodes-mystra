"""Pin the option-to-value tables to the spec and to the published values.

`definition_version` does not hash `decision_values.py`, so editing a table
stales nothing; these tests are what hold it. Run on the laptop with PyYAML:

    python tools/test_decision_values.py
"""
import ast
import re
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src/papers/bmodes/scripts"), str(ROOT / "src/workflow/scripts")]

import catalog_versions
import decision_values as dv
import scale_cuts

ROOT_DECISIONS = yaml.safe_load((ROOT / "astra.yaml").read_text())["decisions"]
TWO_POINT_DECISIONS = yaml.safe_load((ROOT / "analyses/two_point/astra.yaml").read_text())["decisions"]

# The values Paper II was computed with.
PUBLISHED = {
    "fiducial_blind": ("blind_a", "A"),
    "fiducial_catalog": ("size_cut", "SP_v1.4.6.3"),
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
            decisions = TWO_POINT_DECISIONS if name == "patch_count" else ROOT_DECISIONS
            self.assertEqual(set(table), set(decisions[name]["options"]), name)

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
            (ROOT / "analyses/two_point/analyses/xi_scan/xi_scan.sweep.yaml").read_text()
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

    def test_sweep_drivers_refuse_an_empty_version_list(self):
        for driver in ("run_pure_eb_sweep.sh", "run_pure_eb_ptes_sweep.sh"):
            args = ["--versions", "", "--n-samples", "2000", "--integration-method", "m",
                    "--integration-limits", "l", "--rebinning-weights", "w",
                    "--covariance-method", "c"] if driver == "run_pure_eb_sweep.sh" else ["--versions", ""]
            run = subprocess.run(["bash", str(ROOT / "src/papers/bmodes/scripts" / driver), *args],
                                 capture_output=True, text=True)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertIn("--versions is required", run.stderr)


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
