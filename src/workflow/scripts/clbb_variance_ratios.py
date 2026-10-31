"""Measure how much quieter the GLASS mocks are than the data, in C_l^BB.

Two ratios, both against the data's own NaMaster iNKA BB covariance and both
over the harmonic reporting cut 300 <= l <= 1600 (the 17 bandpowers whose
powspace edges fall entirely inside it):

* the **machinery** ratio — the diagonal of the mock-matched BB covariance,
  the paper's own covariance primitive pointed at a mock catalogue, over the
  diagonal of the data's.  This is what the covariance chain predicts the
  mocks' scatter should be, with no rescale hand-fed anywhere;
* the **empirical** ratio — the sample variance of the suite's own BB
  bandpowers over the same data diagonal.  This is what the realisations
  actually do, and it is the truth the machinery ratio is checked against.

The mocks carry no B signal by construction, so both quantities are pure
shape-noise variance ratios and the count mismatch between suite and
catalogue is the whole of what they measure.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plotting_utils import ell_bin_mask, get_powspace_bin_edges  # noqa: E402

ELL_MIN = 300.0
ELL_MAX = 1600.0
BB_ROW = 4  # row 4 of each per-mock Cl block is BB


def load_bb_covariance(path):
    """The BB-BB block of a NaMaster pseudo-Cl bandpower covariance."""

    with fits.open(path) as hdul:
        covariance = np.array(hdul["COVAR_BB_BB"].data, dtype=np.float64)
    if covariance.shape != (32, 32):
        raise ValueError(f"Expected a 32x32 BB covariance in {path}, found {covariance.shape}")
    return covariance


def load_mock_bandpowers(pattern):
    """The suite's BB bandpowers, one row per realisation, plus the l axis."""

    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No per-mock Cl files match {pattern}")
    blocks = [np.load(path) for path in files]
    ell = blocks[0][0]
    return np.array([block[BB_ROW] for block in blocks]), ell, len(files)


def build_ratios(cov_mock, cov_data, bandpowers, ell):
    """Machinery and empirical BB variance ratios on the reporting l cut."""

    mask = ell_bin_mask(ell, ELL_MIN, ELL_MAX)
    if int(mask.sum()) != 17:
        raise ValueError(
            f"Powspace [{ELL_MIN:g}, {ELL_MAX:g}] selected {mask.sum()} bands, expected 17"
        )
    low, high = get_powspace_bin_edges(ell)

    var_data = np.diag(cov_data)[mask]
    var_mock = np.diag(cov_mock)[mask]
    var_empirical = np.var(bandpowers, axis=0, ddof=1)[mask]

    machinery = var_mock / var_data
    empirical = var_empirical / var_data
    return {
        "ell_cut": [ELL_MIN, ELL_MAX],
        "n_bands": int(mask.sum()),
        "n_realisations": int(bandpowers.shape[0]),
        "machinery_variance_ratio_mean": float(np.mean(machinery)),
        "machinery_variance_ratio_median": float(np.median(machinery)),
        "empirical_variance_ratio_mean": float(np.mean(empirical)),
        "empirical_variance_ratio_median": float(np.median(empirical)),
        "per_band": [
            {
                "ell_eff": float(e),
                "ell_low": int(lo),
                "ell_high": int(hi),
                "machinery_variance_ratio": float(m),
                "empirical_variance_ratio": float(p),
            }
            for e, lo, hi, m, p in zip(
                ell[mask], low[mask], high[mask], machinery, empirical
            )
        ],
    }


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cl-cov-mockmatched", required=True)
    parser.add_argument("--cl-cov-data", required=True)
    parser.add_argument("--cl-files", required=True, help="glob for the per-mock Cl blocks")
    parser.add_argument("--out-file", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = _parser(argv)
    bandpowers, ell, n_files = load_mock_bandpowers(args.cl_files)
    result = build_ratios(
        load_bb_covariance(args.cl_cov_mockmatched),
        load_bb_covariance(args.cl_cov_data),
        bandpowers,
        ell,
    )
    result["sources"] = {
        "cl_cov_mockmatched": str(args.cl_cov_mockmatched),
        "cl_cov_data": str(args.cl_cov_data),
        "cl_files": f"{n_files} files matching {args.cl_files}",
    }
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Saved {out_path}: machinery variance ratio "
        f"{result['machinery_variance_ratio_mean']:.4f}, empirical "
        f"{result['empirical_variance_ratio_mean']:.4f} "
        f"({result['n_bands']} bands, {result['n_realisations']} realisations)",
        flush=True,
    )


if __name__ == "__main__":
    main()
