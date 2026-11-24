"""The reference xi_pm covariance on the 20-bin reporting grid.

CosmoCov builds the reference covariance on the 1000-bin integration grid, the
grid every downstream transform wants; the raw xi_pm check wants the 20
reporting bins.

This binds the two: Cov_reporting = B Cov_fine B^T, with B the fine-to-
reporting averaging operator applied to xi_+ and xi_- alike.  TreeCorr's
20-bin xi_pm is a pair-weighted average within each bin, so the default
weighting is the pair density a uniform field gives, w_i ~ theta_r^2 -
theta_l^2 over the fine bins falling inside each reporting bin.  The
unweighted operator -- the ``pointwise`` rebin the pure-mode chain applies to
the *input* xi_pm -- is written alongside as ``cov_pointwise`` so the choice is
visible and cheap to check; the two agree closely, neighbouring fine bins being
almost perfectly correlated.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from decision_values import INTEGRATION_GRID, REPORTING_BINNING


def binning_matrix(left_edges, right_edges, fine_edges, weighting):
    """Fine-to-reporting averaging operator, normalised row by row."""
    n_reporting = len(left_edges)
    n_fine = len(fine_edges) - 1
    fine_centres = np.sqrt(fine_edges[:-1] * fine_edges[1:])

    if weighting == "pairs":
        weights = fine_edges[1:] ** 2 - fine_edges[:-1] ** 2
    elif weighting == "uniform":
        weights = np.ones(n_fine)
    else:
        raise SystemExit(f"unknown weighting {weighting!r}")

    matrix = np.zeros((n_reporting, n_fine))
    for row in range(n_reporting):
        inside = (fine_centres >= left_edges[row]) & (fine_centres < right_edges[row])
        if not inside.any():
            raise SystemExit(
                f"reporting bin {row} ({left_edges[row]:.3f}'-"
                f"{right_edges[row]:.3f}') contains no fine bin; the "
                "covariance grid does not cover the reporting grid."
            )
        matrix[row, inside] = weights[inside] / weights[inside].sum()
    return matrix


def rebin(cov_fine, operator):
    """Apply one operator to both the xi_+ and xi_- halves."""
    n_fine = operator.shape[1]
    block = np.zeros((2 * operator.shape[0], 2 * n_fine))
    block[: operator.shape[0], :n_fine] = operator
    block[operator.shape[0] :, n_fine:] = operator
    return block @ cov_fine @ block.T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--covariance", required=True,
                        help="The fine-grid reference CosmoCov matrix.")
    parser.add_argument("--reporting-binning", required=True,
                        type=REPORTING_BINNING.__getitem__)
    parser.add_argument("--integration-grid", required=True,
                        type=INTEGRATION_GRID.__getitem__)
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args()

    cov_fine = np.loadtxt(args.covariance)
    fine = args.integration_grid
    expected = 2 * fine.nbins
    if cov_fine.shape != (expected, expected):
        raise SystemExit(
            f"{args.covariance} is {cov_fine.shape}, not the "
            f"{expected}x{expected} the declared fine grid implies."
        )

    grid = args.reporting_binning
    edges = np.geomspace(grid.min_sep, grid.max_sep, grid.nbins + 1)
    left, right = edges[:-1], edges[1:]
    fine_edges = np.logspace(np.log10(fine.min_sep),
                             np.log10(fine.max_sep), fine.nbins + 1)

    cov = rebin(cov_fine, binning_matrix(left, right, fine_edges, "pairs"))
    cov_pointwise = rebin(cov_fine,
                          binning_matrix(left, right, fine_edges, "uniform"))

    nbins = len(left)
    sigma = np.sqrt(np.diag(cov))
    np.savez(
        args.out_file,
        cov=cov,
        cov_pointwise=cov_pointwise,
        sigma_xip=sigma[:nbins],
        sigma_xim=sigma[nbins:],
        left_edges=left,
        right_edges=right,
        provenance=json.dumps({
            "covariance": args.covariance,
            "fine_grid": fine._asdict(),
            "weighting": "pairs",
        }),
    )
    ratio = np.sqrt(np.diag(cov_pointwise)) / sigma
    print(f"Saved {args.out_file}")
    print(f"pointwise/pair-weighted sigma ratio: "
          f"{ratio.min():.4f}-{ratio.max():.4f}")
    for i in range(nbins):
        print(f"{0.5 * (left[i] + right[i]):8.2f}  sigma_xip {sigma[i]:11.4e}  "
              f"sigma_xim {sigma[nbins + i]:11.4e}")


if __name__ == "__main__":
    main()
