"""A GLASS suite's density and per-component σ_e, as JSON.

    python glass_mock_sampling.py --suite SUITE --fiducial-catalog CAT --summary SUMMARY --nz NZ --out-file S.json
"""

import argparse

from glass_mock_generate import sampling

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--suite", required=True)
p.add_argument("--fiducial-catalog", required=True)
p.add_argument("--summary", required=True)
p.add_argument("--nz", required=True)
p.add_argument("--out-file", required=True)
a = p.parse_args()
sampling(a.suite, a.fiducial_catalog, a.summary, a.nz, a.out_file)
