"""Bundle a reference catalog's fiducial-chain artifacts for the comparison readers.

A catalog comparison sets the universe's fiducial catalog against a reference
catalog whose products were made as *another* universe's fiducial chain (the
published v1.4.6.3, materialized in `baseline`). The version-comparison and
PTE-matrix readers find every non-fiducial version by its July-canonical file
name inside one directory (`--results-dir`, `--cov-dir`, `--pte-intermediate-dir`,
`--cosebis-pte-dir`), so this writes all of the reference's artifacts under those
names into one tar.gz that every reader can untar and point at.

The reference version is the config's single non-fiducial version.
"""

import argparse
import os
import shutil
import tarfile
import tempfile

from bundle_sweep import BLIND, BINNING, COV_GRID, XI_GRIDS, XI_NPATCH
from paper_config import load_config
from sweep_versions import nonfiducial_versions


def canonical_names(ver):
    """Argument name → the July-canonical file name its artifact is read under."""
    xi = {
        f"xi_{grid}": (
            f"{ver}_xi_minsep={g['min_sep']}_maxsep={g['max_sep']}"
            f"_nbins={g['nbins']}_npatch={XI_NPATCH}.txt"
        )
        for grid, g in XI_GRIDS.items()
    }
    cov = (
        f"covariance_{ver}_{BLIND}_g_minsep={COV_GRID['min_sep']}"
        f"_maxsep={COV_GRID['max_sep']}_nbins={COV_GRID['nbins']}_masked"
    )
    cl = {
        f"{arg}": f"pseudo_{product}_{ver}_blind={BLIND}_{BINNING}_nbins={nbins}.fits"
        for arg, product, nbins in [
            ("cl", "cl", 32),
            ("cl_cov", "cl_cov", 32),
            ("cl_fine", "cl", 96),
            ("cl_cov_fine", "cl_cov", 96),
        ]
    }
    return {
        "xi_reporting": xi["xi_reporting"],
        "xi_integration": xi["xi_integration"],
        "cov_integration": os.path.join(cov, f"{cov}_processed.txt"),
        **cl,
        "pure_eb": f"{ver}_{BLIND}_pure_eb_semianalytic.npz",
        "pure_eb_ptes": f"{ver}_{BLIND}_pure_eb_ptes.npz",
        "cosebis_ptes": f"cosebis_ptes_{ver}_{BLIND}.npz",
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, help="Catalog config (config/catalogs/<option>.yaml)")
    for arg in canonical_names("V"):
        ap.add_argument(f"--{arg.replace('_', '-')}", dest=arg, required=True)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    (ver,) = nonfiducial_versions(load_config(a.config))
    scratch = tempfile.mkdtemp(prefix="bundle_reference_", dir=os.environ.get("TMPDIR"))
    try:
        for arg, name in canonical_names(ver).items():
            dest = os.path.join(scratch, name)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copyfile(getattr(a, arg), dest)
            print(f"[bundle_reference] {getattr(a, arg)} -> {name}")
        with tarfile.open(a.out_file, "w:gz") as tar:
            tar.add(scratch, arcname=".")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    print(f"[bundle_reference] wrote {a.out_file} ({ver})")


if __name__ == "__main__":
    main()
