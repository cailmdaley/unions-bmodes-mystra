"""CAMB matter Cℓ of every pair of GLASS shells, in the suite's cosmology.

    python glass_mock_shell_cls.py --out-file CLS.npy
"""

import argparse

from glass_mock_generate import shell_cls

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--out-file", required=True)
a = p.parse_args()
shell_cls(a.out_file)
