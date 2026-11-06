"""Calculate PTE matrices for pure E/B-mode scale-cut robustness.

Reads the gathered
pure-E/B ``semianalytic.npz`` (data vectors + MC covariance), evaluates the
ξ_+^B / ξ_-^B / joint ξ_tot^B χ² PTE matrices over the scale-cut grid via
``sp_validation.b_modes.calculate_eb_statistics``, and writes the PTE matrices either to a caller-given path
(``--out-file``) or, in directory mode (``--out``, used by the sweep driver),
to ``{out}/{version}_{blind}_pure_eb_ptes.npz``.

The Hartlap factor follows the covariance artifact
(``pure_mode_covariance.hartlap_samples``): the Monte Carlo draw count it
records, or no correction for an exact covariance.
Artifacts that record no draw count use the shared 2000-draw convention.

    python calculate_pure_eb_ptes.py \
        --version SP_v1.4.6.3_leak_corr --blind A \
        --pure-eb-data <..._pure_eb_semianalytic.npz> \
        --cov-integration <cov ..._processed.txt> \
        --npatch 1 --out-file <ptes.npz>
"""

import argparse
import os

import numpy as np
from pure_mode_covariance import hartlap_samples

from sp_validation.b_modes import calculate_eb_statistics

#: calculate_eb_statistics always multiplies by (N − p − 2)/(N − 1); at this N
#: that ratio rounds to exactly 1.0, which is how "no correction" reaches it.
NO_HARTLAP_N = 10**300


class FakeGG:
    """Minimal GGCorrelation-like object for calculate_eb_statistics."""

    def __init__(self, nbins, npatch):
        self.nbins = nbins
        self.npatch1 = npatch
        self.npatch2 = npatch


def calculate_ptes(
    version,
    blind,
    pure_eb_data,
    cov_integration,
    npatch,
    output_dir=None,
    out_file=None,
):
    dataset = np.load(pure_eb_data)

    theta = dataset["theta"]
    nbins = len(theta)

    results = {
        "gg": FakeGG(nbins, int(npatch)),
        "xip_E": dataset["xip_E"],
        "xim_E": dataset["xim_E"],
        "xip_B": dataset["xip_B"],
        "xim_B": dataset["xim_B"],
        "xip_amb": dataset["xip_amb"],
        "xim_amb": dataset["xim_amb"],
        "cov": dataset["cov_pure_eb"],
    }

    n_hartlap = hartlap_samples(dataset)
    print(f"Calculating PTE matrices for {version} (Hartlap N={n_hartlap})...")
    results = calculate_eb_statistics(
        results,
        cov_path_int=cov_integration,
        n_samples=NO_HARTLAP_N if n_hartlap is None else n_hartlap,
    )

    pte_matrices = results["pte_matrices"]
    output_data = {
        "theta": theta,
        "pte_xip_B": pte_matrices["xip_B"],
        "pte_xim_B": pte_matrices["xim_B"],
        "pte_combined": pte_matrices["combined"],
    }

    if out_file is not None:
        out_path = out_file
    else:
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, f"{version}_{blind}_pure_eb_ptes.npz")
    np.savez(out_path, **output_data)
    print(f"Saved PTE matrices to {out_path}")
    return out_path


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--version", required=True)
    ap.add_argument("--blind", default="A")
    ap.add_argument("--pure-eb-data", required=True, help="Gathered semianalytic .npz")
    ap.add_argument(
        "--cov-integration",
        default=None,
        help="Integration-grid covariance _processed.txt (optional)",
    )
    ap.add_argument("--npatch", type=int, default=1)
    ap.add_argument(
        "--out",
        default=None,
        help="Output directory; writes {version}_{blind}_pure_eb_ptes.npz "
        "(sweep-driver mode)",
    )
    ap.add_argument(
        "--out-file", default=None, help="Output NPZ path (lc {output})"
    )
    a = ap.parse_args(argv)
    if (a.out is None) == (a.out_file is None):
        ap.error("exactly one of --out / --out-file is required")
    calculate_ptes(
        version=a.version,
        blind=a.blind,
        pure_eb_data=a.pure_eb_data,
        cov_integration=a.cov_integration,
        npatch=a.npatch,
        output_dir=a.out,
        out_file=a.out_file,
    )


if __name__ == "__main__":
    _from_cli()
