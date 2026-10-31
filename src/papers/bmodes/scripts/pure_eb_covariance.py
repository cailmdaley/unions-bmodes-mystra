"""Pure E/B covariance structure claim.

Validates covariance matrix for B-mode tests by analyzing block structure:
- 6 blocks: E+, E-, B+, B-, amb+, amb-
- Block-wise condition numbers
- Correlation matrix heatmap

Shows that E and B blocks are well-conditioned while ill-conditioning
is localized to ambiguous modes (expected from finite angular range).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from plotting_utils import FIG_WIDTH_SINGLE, PAPER_MPLSTYLE
from paper_config import load_config

plt.style.use(PAPER_MPLSTYLE)


def _analyze_block(covariance, start_idx, end_idx, name):
    """Compute condition number and definiteness for a covariance block.

    Parameters
    ----------
    covariance : array
        Full covariance matrix
    start_idx : int
        Start index of block
    end_idx : int
        End index of block (exclusive)
    name : str
        Block name for reporting

    Returns
    -------
    dict
        Block statistics: condition number, positive definiteness, eigenvalue bounds
    """
    block = covariance[start_idx:end_idx, start_idx:end_idx]
    eigenvalues = np.linalg.eigvalsh(block)

    min_eig = float(eigenvalues[0])
    max_eig = float(eigenvalues[-1])
    is_positive_definite = bool(min_eig > 0)
    condition_number = max_eig / min_eig if min_eig > 0 else float("inf")

    return {
        "condition_number": condition_number,
        "is_positive_definite": is_positive_definite,
        "min_eigenvalue": min_eig,
        "max_eigenvalue": max_eig,
        "block_size": list(block.shape),
        "description": name,
    }


def _cov_to_corr(covariance):
    """Convert covariance matrix to correlation matrix.

    Parameters
    ----------
    covariance : array (N, N)
        Covariance matrix

    Returns
    -------
    array (N, N)
        Correlation matrix with diagonal = 1
    """
    std = np.sqrt(np.diag(covariance))
    std_inv = 1.0 / np.where(std > 0, std, 1.0)
    correlation = covariance * np.outer(std_inv, std_inv)
    return correlation


def main(config, pure_eb_path, out_file):
    version = config["fiducial"]["version"]

    # Load precomputed pure E/B data
    dataset = np.load(pure_eb_path)

    theta = dataset["theta"]
    nbins = len(theta)
    cov_pure_eb = dataset["cov_pure_eb"]

    # Matrix is 6*nbins x 6*nbins:
    # Block structure: [E+, E-, B+, B-, amb+, amb-]
    # Each block is nbins x nbins
    assert cov_pure_eb.shape == (6 * nbins, 6 * nbins), (
        f"Expected ({6 * nbins}, {6 * nbins}), got {cov_pure_eb.shape}"
    )

    # Analyze full matrix
    eigenvalues = np.linalg.eigvalsh(cov_pure_eb)
    min_eig = float(eigenvalues[0])
    max_eig = float(eigenvalues[-1])
    is_positive_definite = bool(min_eig > 0)
    condition_number = max_eig / min_eig if min_eig > 0 else float("inf")

    # Analyze blocks
    block_analysis = {
        "xi_E": _analyze_block(
            cov_pure_eb, 0, 2 * nbins, "xi_+^E and xi_-^E combined (40x40)"
        ),
        "xi_B": _analyze_block(
            cov_pure_eb, 2 * nbins, 4 * nbins, "xi_+^B and xi_-^B combined (40x40)"
        ),
        "xi_amb": _analyze_block(
            cov_pure_eb, 4 * nbins, 6 * nbins, "xi_+^amb and xi_-^amb combined (40x40)"
        ),
    }

    # Convert to correlation matrix for visualization
    correlation = _cov_to_corr(cov_pure_eb)

    # Create figure (matching v1 epistemics style)
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_SINGLE, FIG_WIDTH_SINGLE * 0.8))

    # Plot correlation matrix (no origin="lower" so y increases downward like v1)
    im = ax.imshow(
        correlation,
        cmap="icefire",
        vmin=-1.0,
        vmax=1.0,
        aspect="equal",
    )

    # 6 block labels (matching v1 epistemics style)
    block_labels = [
        r"$\xi_+^{\mathrm{E}}$",
        r"$\xi_-^{\mathrm{E}}$",
        r"$\xi_+^{\mathrm{B}}$",
        r"$\xi_-^{\mathrm{B}}$",
        r"$\xi_+^{\mathrm{amb}}$",
        r"$\xi_-^{\mathrm{amb}}$",
    ]

    # Compute block boundaries and centers
    block_sizes = [nbins] * 6
    boundaries = np.cumsum(block_sizes)
    centers = [
        0.5 * (start + end) - 0.5
        for start, end in zip([0] + list(boundaries[:-1]), boundaries)
    ]

    # Set ticks at block centers with labels (no axis labels, tick labels only)
    ax.set_xticks(centers)
    ax.set_xticklabels(block_labels)
    ax.set_yticks(centers)
    ax.set_yticklabels(block_labels, rotation=0)

    # Colorbar — match height to matrix
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Correlation coefficient")

    fig.tight_layout()

    # Save output: a single 300-dpi PNG at the caller-given path
    out_file = Path(out_file)
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"Saved {out_file}")

    plt.close(fig)

    # Print summary
    print("\nCovariance Matrix Summary:")
    print(f"  Shape: {cov_pure_eb.shape}")
    print(f"  Full matrix condition number: {condition_number:.2e}")
    print(f"  Positive definite: {is_positive_definite}")
    print("\nBlock Analysis:")
    for block_name, stats in block_analysis.items():
        print(f"  {block_name}:")
        print(f"    Condition number: {stats['condition_number']:.2e}")
        print(f"    Positive definite: {stats['is_positive_definite']}")
        print(f"    Min eigenvalue: {stats['min_eigenvalue']:.2e}")
        print(f"    Max eigenvalue: {stats['max_eigenvalue']:.2e}")


def _from_cli(argv=None):
    import yaml

    ap = argparse.ArgumentParser(
        description="Pure E/B 6-block covariance correlation-matrix paper figure (fiducial catalog)."
    )
    ap.add_argument(
        "--config", required=True, help="Absolute path to bmodes config.yaml"
    )
    ap.add_argument(
        "--pure-eb-data",
        required=True,
        help="Fiducial <version>_<blind>_pure_eb_semianalytic.npz "
        "(provides the 6-block cov_pure_eb)",
    )
    ap.add_argument(
        "--out-file", required=True, help="Output PNG path (lc {output})"
    )
    a = ap.parse_args(argv)
    config = load_config(a.config)
    main(config, a.pure_eb_data, a.out_file)


if __name__ == "__main__":
    _from_cli()
