"""Parameter-plumbing regression tests; run in the project environment on Slurm.

These use small synthetic arrays and process stubs, not the GLASS ensemble.
Without `sp_validation` the module skips.
"""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

if importlib.util.find_spec("sp_validation") is None:
    raise unittest.SkipTest("needs the project environment")

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src/papers/bmodes/scripts"), str(ROOT / "src/workflow/scripts")]

import bundle_sweep
from catalog_versions import VERSIONS, Catalogs
import compute_cosebis_pte_single as pairs
from decision_values import MC_SAMPLES
import generate_cl as cl
import generate_cl_cov as cl_cov
import merge_cosebis_pte_rows as merge
import primary_null_test_report as reporting


class DecisionWiringTests(unittest.TestCase):
    def test_scale_cut_pair_options_and_row_partition(self):
        production = pairs._pte_scale_cut_pairs()
        self.assertEqual(len(production), 206)
        self.assertNotIn((9, 10), production)
        for expected in [production]:
            rows = [pairs._pte_scale_cut_pairs(row=i) for i in range(20)]
            self.assertEqual([pair for row in rows for pair in row], expected)
            with tempfile.TemporaryDirectory() as tmp:
                members = []
                for i, row in enumerate(rows):
                    p = Path(tmp) / f"row_{i}.npz"
                    arrays = {k: np.zeros((21, 21)) for k in merge.MATRICES}
                    np.savez(p, version="test", blind="A", nmodes=20,
                             mode_subsets=np.array([6, 20]), theta_grid=np.arange(21),
                             pairs=np.asarray(row, dtype=int).reshape(-1, 2),
                             theta_min=np.array([x for x, y in row]),
                             theta_max=np.array([y for x, y in row]),
                             En=np.zeros((len(row), 20)), Bn=np.zeros((len(row), 20)),
                             convergence_failure=np.zeros(len(row), dtype=bool), **arrays)
                    members.append(str(p))
                output = Path(tmp) / "merged.npz"
                merge.main(["--members", *members, "--out-file", str(output)])
                with np.load(output) as data:
                    np.testing.assert_array_equal(data["pairs"], expected)
                    self.assertEqual(data["En"].shape, (len(expected), 20))

    def test_harmonic_count_changes_bundle_names_not_member_ids(self):
        versions = ["SP_v1.4.5", "SP_v1.4.5_leak_corr"]
        catalogs = Catalogs("size_cut", "applied")
        default = bundle_sweep._cl_map(versions, 96, catalogs)
        alternative = bundle_sweep._cl_map(versions, 128, catalogs)
        self.assertEqual(set(default), set(alternative))
        for key in default:
            self.assertEqual(alternative[key], default[key].replace("nbins=96", "nbins=128"))

    def test_estimator_and_noise_reach_existing_primitives(self):
        for method in ["catalog"]:
            with tempfile.TemporaryDirectory() as tmp, patch.object(cl, "CosmologyValidation") as cv:
                cl.generate_cl("test", tmp, "unused", binning="powspace", nbins=128)
                self.assertEqual(cv.call_args.kwargs["cell_method"], method)
                self.assertEqual(cv.call_args.kwargs["n_ell_bins"], 128)
                cv.return_value.calculate_pseudo_cl.assert_called_once()
            with tempfile.TemporaryDirectory() as tmp, patch.object(cl_cov, "CosmologyValidation") as cv:
                cl_cov.generate_cl_cov("test", tmp, "unused", nbins=128)
                cv.return_value.calculate_pseudo_cl_eb_cov.assert_called_once()

    def test_bin_count_reaches_the_measurement(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = str(Path(tmp) / "cl.fits")
            def fake_generate(**kwargs):
                self.assertEqual(kwargs["nbins"], 128)
                path = Path(kwargs["output_dir"]) / "generated.fits"
                path.write_bytes(b"fixture")
                return str(path)
            with patch.object(cl, "generate_cl", side_effect=fake_generate):
                cl._from_cli(["--ver", "test", "--cat-config", "unused", "--nbins", "128",
                              "--out-file", output])
            self.assertEqual(Path(output).read_bytes(), b"fixture")

    def test_leakage_selects_primary_and_complete_scan_twins(self):
        original = Catalogs("size_cut", "applied")
        raw = Catalogs("size_cut", "uncorrected")
        self.assertEqual(raw.fiducial, "SP_v1.4.6.3")
        self.assertEqual(raw.comparison, [v.removesuffix("_leak_corr") for v in original.comparison])
        self.assertEqual(raw.labels, {k.removesuffix("_leak_corr"): v for k, v in original.labels.items()})
        slots = original.nonfiducial
        self.assertEqual(raw.scan_slots, slots)
        for catalogs in (original, raw):
            actual = [catalogs.scan_version(v) for v in slots]
            self.assertEqual(set(actual + [catalogs.fiducial]), {v for v in VERSIONS if "_ecut" not in v})
            for mapping in (bundle_sweep._xi_map(slots, catalogs), bundle_sweep._cov_map(slots, catalogs),
                            bundle_sweep._cosebis_ptes_map(slots, catalogs), bundle_sweep._cl_map(slots, 128, catalogs)):
                self.assertEqual(len(set(mapping.values())), len(mapping))
            fine = bundle_sweep._cl_map(slots, 128, catalogs)
            for name, path in fine.items():
                if name.endswith("_fine"):
                    self.assertEqual("_leak_corr" in path, catalogs is original)
        self.assertEqual(raw.scan_version("SP_v1.4.6.3"), "SP_v1.4.6.3_leak_corr")

    def test_selected_family_is_recalibrated_from_selected_columns(self):
        names = ["B_n<=6", "B_n<=20", "xi+^B", "xi-^B", "xi_tot^B", "C_ell^BB"]
        ptes = np.tile([.4, .001, .5, .6, .7, .8], (10, 1))
        reference = dict(observed_ptes=dict(zip(names, [.3, .002, .4, .5, .6, .7])),
                         window_bins=[9, 15], scale_cut_arcmin=[12, 83],
                         scale_cut_convention="grid_bin_edges", global_alpha=.05)
        full = dict(reference, window_bins=[0, 19], scale_cut_arcmin=[1, 250])
        reference["full_range"] = full
        both = reporting.select_report(reference, names, ptes, "adopted", "n_le_6_and_n_le_20")
        six = reporting.select_report(reference, names, ptes, "full_range", "n_le_6_only")
        self.assertAlmostEqual(both["calibration"]["threshold_one_sided"], .001)
        self.assertAlmostEqual(six["calibration"]["threshold_one_sided"], .4)
        self.assertEqual(six["calibration"]["window_bins"], [0, 19])
        self.assertNotIn("B_n<=20", six["calibration"]["stat_names"])

    def test_harmonic_alternative_keeps_six_modes(self):
        import harmonic_config_cosebis_comparison as harmonic
        self.assertEqual(harmonic._RELIABLE_MODE_MAX_BY_NBINS[128], 6)

    def test_mc_option_reaches_all_chunks(self):
        import calculate_pure_eb_ptes as ptes
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            artifact = tmp / "covariance.npz"
            vectors = {k: np.zeros(2) for k in
                       ("theta", "xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")}
            for method, count in [("monte_carlo", 500), ("monte_carlo", 1000),
                                  ("monte_carlo", 2000), ("exact_linear", None)]:
                record = dict(covariance_method=method)
                if count is not None:
                    record["n_samples"] = count
                np.savez(artifact, cov_pure_eb=np.eye(12), **vectors, **record)
                matrices = dict(xip_B=np.zeros((2, 2)), xim_B=np.zeros((2, 2)), combined=np.zeros((2, 2)))
                with patch.object(ptes, "calculate_eb_statistics", return_value=dict(pte_matrices=matrices)) as stats:
                    ptes.calculate_ptes("test", "A", artifact, artifact, 1, out_file=tmp / "ptes.npz")
                    self.assertEqual(stats.call_args.kwargs["n_samples"],
                                     ptes.NO_HARTLAP_N if count is None else count)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            python = tmp / "python"
            python.write_text(f"#!{sys.executable}\n" + '''import json, pathlib, sys
args = sys.argv[2:]
def get(flag): return args[args.index(flag)+1]
if pathlib.Path(sys.argv[1]).name == "precompute_pure_eb_chunk.py":
    out = pathlib.Path(get("--out")); out.mkdir(parents=True, exist_ok=True)
    (out / ("pure_eb_chunk_"+get("--chunk-id")+".npz")).write_text(get("--n-samples"))
else:
    chunks = pathlib.Path(get("--chunks-dir"))
    values = [int(p.read_text()) for p in sorted(chunks.glob("*.npz"))]
    pathlib.Path(get("--out-file")).write_text(json.dumps(values))
''')
            python.chmod(0o755)
            env = dict(os.environ, PATH=f"{tmp}:{os.environ['PATH']}")
            for option, expected in MC_SAMPLES.items():
                output = tmp / f"{option}.npz"
                subprocess.run(["bash", str(ROOT / "src/papers/bmodes/scripts/run_pure_eb_covariance.sh"),
                                "--version", "test", "--blind", "A", "--cat-config", "unused",
                                "--xi-reporting", "unused", "--xi-integration", "unused",
                                "--cov-integration", "unused", "--integration-method", "adaptive_dqags",
                                "--integration-limits", "reporting_extent", "--rebinning-weights", "none",
                                "--covariance-method", "monte_carlo", "--n-samples", str(expected),
                                "--n-chunks", "2", "--nproc", "2", "--out-file", str(output)],
                               check=True, env=env)
                self.assertEqual(json.loads(output.read_text()), [expected, expected])

    def test_cosmocov_options_reach_real_ini_generator(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cosmology = tmp / "cosmology.json"
            cosmology.write_text(json.dumps(dict(Omega_m=.3, Omega_v=.7, sigma_8=.8,
                                                 n_s=.96, Omega_b=.05, h=.7)))
            mask = tmp / "mask.txt"
            mask.write_text("0 1\n")
            binary = tmp / "cosmocov"
            binary.write_text(f"#!{sys.executable}\n" + '''import pathlib, sys
ini = dict(line.split(":", 1) for line in pathlib.Path(sys.argv[2]).read_text().splitlines() if ":" in line)
ini = {k.strip(): v.strip() for k, v in ini.items()}
i = int(sys.argv[1]); pm = ["++", "--", "+-"][i-1]
pathlib.Path(f"cov_tmp_ssss_{pm}_cov_Ntheta{ini['ntheta']}_Ntomo1_{i}").write_text(f"{i},{ini['ng']},{ini['c_footprint_file']}\\n")
''')
            binary.chmod(0o755)
            for cls, ng in [("gaussian_only", "0"), ("ng", "1")]:
                for treatment in [None, "masked", "unmasked"]:
                    output = tmp / f"{cls}_{treatment}.txt"
                    subprocess.run(["bash", str(ROOT / "src/workflow/scripts/run_cosmocov_chain.sh"),
                                    "--version", "SP_v1.4.6.3", "--blind", "blind_a", "--min-sep", "1",
                                    "--max-sep", "250", "--nbins", "3", "--gaussian", cls,
                                    *(["--mask-application", treatment] if treatment else []),
                                    "--mask-cls", str(mask),
                                    "--planck18-json", str(cosmology),
                                    "--cosmocov", str(binary), "--libgsl", str(mask),
                                    "--emit", "blocks", "--out-file", str(output)], check=True)
                    footprint = "mask_cls.txt" if treatment != "unmasked" and ng == "1" else "no_footprint_mask"
                    self.assertEqual(output.read_text(), "".join(f"{i},{ng},{footprint}\n" for i in (1, 2, 3)))


if __name__ == "__main__":
    unittest.main()
