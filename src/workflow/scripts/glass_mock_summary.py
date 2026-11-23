"""The catalog summary of one GLASS realization, its catalog discarded.

Generates realization N in a scratch directory (``glass_mock_generate.py``) and
writes its ``catalog_summary.py`` entry under KEY, the cat_config entry whose
n(z) the Cℓ covariance of that catalog takes.

    python glass_mock_summary.py --number N --key KEY --summary SUMMARY --suite SUITE \
        --fiducial-catalog CAT --nz NZ --shell-cls C --sampling S --out-file realization.npz
"""

import argparse
import os
import tempfile
from pathlib import Path

import numpy as np

from catalog_summary import summarize
from glass_mock_generate import footprint_mask, realization


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--number", type=int, required=True)
    p.add_argument("--key", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--suite", required=True)
    p.add_argument("--fiducial-catalog", required=True)
    p.add_argument("--nz", required=True)
    p.add_argument("--shell-cls", required=True)
    p.add_argument("--sampling", required=True)
    p.add_argument("--out-file", required=True)
    a = p.parse_args(argv)
    with tempfile.TemporaryDirectory(dir=os.environ.get("SLURM_TMPDIR")) as tmp:
        catalog = Path(tmp) / f"unions_glass_sim_{a.number:05d}_4096.fits"
        mask = footprint_mask(a.summary, a.suite, a.fiducial_catalog)
        realization(a.number, mask, a.nz, a.shell_cls, a.sampling, catalog)
        entry = summarize(str(catalog))
    with open(a.out_file, "wb") as f:
        np.savez_compressed(f, **{f"{a.key}/{name}": value for name, value in entry.items()})


if __name__ == "__main__":
    main()
