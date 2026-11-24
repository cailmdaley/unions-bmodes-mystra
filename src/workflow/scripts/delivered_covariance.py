"""One delivered CosmoCov ξ± covariance, as the text matrix the analysis reads.

The `xi_covariances` input is a tar holding, per catalog version,
`integration/<version>.npy` (the upper triangle of the Gaussian 1000-bin matrix,
row by row) and `reporting_ng/<version>.txt` (the Gaussian plus non-Gaussian
20-bin matrix). A catalog's corrected and uncorrected forms share one entry.

    python delivered_covariance.py --covariances <xi_covariances.tar> \\
        --grid integration --fiducial-catalog initial --out-file <cov.txt>
"""

import argparse
import io
import tarfile

import numpy as np
from decision_values import FIDUCIAL_CATALOG


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--covariances", required=True)
    ap.add_argument("--grid", choices=["integration", "reporting_ng"], required=True)
    ap.add_argument("--fiducial-catalog", type=FIDUCIAL_CATALOG.__getitem__, required=True)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)
    suffix = {"integration": "npy", "reporting_ng": "txt"}[a.grid]
    with tarfile.open(a.covariances) as tar:
        data = tar.extractfile(f"{a.grid}/{a.fiducial_catalog}.{suffix}").read()
    if a.grid == "reporting_ng":
        with open(a.out_file, "wb") as f:
            f.write(data)
        return
    triangle = np.load(io.BytesIO(data))
    n = int((np.sqrt(8 * triangle.size + 1) - 1) / 2)
    cov = np.zeros((n, n))
    cov[np.triu_indices(n)] = triangle
    np.savetxt(a.out_file, cov + np.triu(cov, 1).T)


if __name__ == "__main__":
    _from_cli()
