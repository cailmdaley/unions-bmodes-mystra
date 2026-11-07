"""Split files recombine into the recorded npz, and hold what each side should.

    python tools/test_split_covariances.py
"""
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tools"), str(ROOT / "src/papers/bmodes/scripts")]

import pure_mode_covariance
import split_recorded_covariances as splitter

VECTORS = ("xip_total", "xim_total", "xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")


def recorded_pure_eb(rng, **record):
    arrays = {k: rng.normal(size=20) for k in VECTORS}
    return {"theta": rng.normal(size=20), "theta_int": rng.normal(size=1000),
            **arrays, "cov_pure_eb": rng.normal(size=(120, 120)), **record}


def recorded_cosebis(rng):
    arrays = {"nmodes": 20, "version": "SP_v1.4.6.3_leak_corr"}
    for cut in ("fiducial", "full"):
        arrays |= {f"{cut}_En": rng.normal(size=20), f"{cut}_Bn": rng.normal(size=20),
                   f"{cut}_cov": rng.normal(size=(40, 40)), f"{cut}_sigma_E": rng.normal(size=20),
                   f"{cut}_sigma_B": rng.normal(size=20), f"{cut}_scale_cut": rng.normal(size=2)}
    return arrays


class Split(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(1)
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def split(self, kind, arrays):
        recorded = self.dir / "recorded.npz"
        np.savez(recorded, **arrays)
        data, covariance = self.dir / "data.npz", self.dir / "cov.npz"
        splitter.split_file(kind, recorded, data, covariance)
        return np.load(data), np.load(covariance)

    def assert_recombines(self, arrays, data, covariance):
        merged = splitter.reassemble(dict(data), dict(covariance))
        self.assertEqual(set(merged), set(arrays))
        for key, value in arrays.items():
            np.testing.assert_array_equal(merged[key], value, key)

    def test_pure_eb_exact(self):
        arrays = recorded_pure_eb(self.rng, covariance_method="exact_linear")
        data, covariance = self.split("pure_eb", arrays)
        self.assertEqual(set(data.files), {"theta", "theta_int", *VECTORS})
        self.assertEqual(set(covariance.files), {"theta", "theta_int", "cov_pure_eb", "covariance_method"})
        self.assertIsNone(pure_mode_covariance.hartlap_samples(covariance))
        self.assert_recombines(arrays, data, covariance)

    def test_pure_eb_monte_carlo_keeps_its_hartlap_record(self):
        arrays = recorded_pure_eb(self.rng, covariance_method="monte_carlo", n_samples=2000)
        data, covariance = self.split("pure_eb", arrays)
        self.assertEqual(pure_mode_covariance.hartlap_samples(covariance), 2000)
        self.assert_recombines(arrays, data, covariance)

    def test_cosebis(self):
        arrays = recorded_cosebis(self.rng)
        data, covariance = self.split("cosebis", arrays)
        self.assertFalse({k for k in data.files if "cov" in k or "sigma" in k})
        self.assertFalse({k for k in covariance.files if k.endswith(("En", "Bn"))})
        self.assert_recombines(arrays, data, covariance)

    def test_pure_eb_sweep_names_each_side(self):
        member = self.dir / "SP_v1.4.5_A_pure_eb_semianalytic.npz"
        np.savez(member, **recorded_pure_eb(self.rng, covariance_method="exact_linear"))
        bundle = self.dir / "sweep.tar.gz"
        with tarfile.open(bundle, "w:gz") as tar:
            tar.add(member, arcname=member.name)
        splitter.split_sweep(bundle, self.dir / "d.tar.gz", self.dir / "c.tar.gz")
        with tarfile.open(self.dir / "d.tar.gz") as data, tarfile.open(self.dir / "c.tar.gz") as cov:
            self.assertIn("./SP_v1.4.5_pure_eb_data.npz", data.getnames())
            self.assertIn("./SP_v1.4.5_A_pure_eb_covariance.npz", cov.getnames())


if __name__ == "__main__":
    unittest.main()
