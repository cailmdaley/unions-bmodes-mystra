"""Gather the GLASS suite members into the three ensembles.

Each member is a ``glass_mock_member.py`` tar.gz. Writes one tar.gz holding
``mock_xi_fine.npz`` and ``mock_xi_reporting.npz`` (``mock_xi_ensemble.stack``
layout: ``(n_mocks, nbins)`` blocks of meanr, xip, xim, npairs, weight, with
the mock IDs and bin edges) and ``cl/cl_glass_mock_<id>_4096.npy``.

    python glass_mock_suite.py --members m1.tar.gz ... --out-file suite.tar.gz
    python glass_mock_suite.py --extract mock_xi_fine.npz --suite suite.tar.gz --out-file x.npz
"""

import argparse
import glob
import shutil
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from decision_values import INTEGRATION_GRID, REPORTING_BINNING
from mock_xi_ensemble import stack


def gather(members, out_file):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for member in members:
            with tarfile.open(member) as tar:
                tar.extractall(tmp, filter="data")
        for sub, grid, name in (
            ("xi_fine", INTEGRATION_GRID["bins_1000"], "mock_xi_fine.npz"),
            ("xi_reporting", REPORTING_BINNING["bins_20"], "mock_xi_reporting.npz"),
        ):
            files = sorted(glob.glob(str(tmp / sub / "*.fits")))
            np.savez(tmp / name, **stack(files, min_sep=grid.min_sep, max_sep=grid.max_sep, nbins=grid.nbins))
        n_cl = len(list((tmp / "cl").glob("*.npy")))
        if not len(members) == n_cl:
            raise SystemExit(f"{len(members)} members but {n_cl} Cℓ files")
        with tarfile.open(out_file, "w:gz") as tar:
            for name in ("mock_xi_fine.npz", "mock_xi_reporting.npz", "cl"):
                tar.add(tmp / name, arcname=name)


def extract(suite, member, out_file):
    """Copy one entry of the suite out: an npz, or the cl/ directory as a tar.gz."""
    with tarfile.open(suite) as tar, tempfile.TemporaryDirectory() as tmp:
        if member == "cl":
            tar.extractall(tmp, members=[m for m in tar.getmembers() if m.name.startswith("cl")], filter="data")
            with tarfile.open(out_file, "w:gz") as out:
                for path in sorted(Path(tmp, "cl").glob("*.npy")):
                    out.add(path, arcname=path.name)
        else:
            tar.extract(member, tmp, filter="data")
            shutil.move(Path(tmp, member), out_file)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--members", nargs="+")
    p.add_argument("--extract", choices=["mock_xi_fine.npz", "mock_xi_reporting.npz", "cl"])
    p.add_argument("--suite")
    p.add_argument("--out-file", required=True)
    a = p.parse_args(argv)
    if a.extract:
        extract(a.suite, a.extract, a.out_file)
    else:
        gather(a.members, a.out_file)


if __name__ == "__main__":
    main()
