"""Figure for the C_l^BB covariance E-to-B mixing-term test.

Upper panel: sigma(C_l^BB) per bandpower for the suite's empirical scatter and
for the three analytic variants (chain input, E-modes zeroed, mock mean E-mode
input). Lower panel: sigma_emp / sigma_ana for each variant, with unity marked.
The harmonic reporting cut 300 < l < 1600 is shaded in both.
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from plotting_utils import (  # noqa: E402
    CATEGORICAL_COLORS,
    FIG_WIDTH_FULL,
    PAPER_MPLSTYLE,
)

VARIANTS = [
    ("a_theory", r"analytic, chain input $C_\ell^{EE}$", CATEGORICAL_COLORS[0], "-"),
    ("b_noEE", r"analytic, $C_\ell^{EE} = 0$", CATEGORICAL_COLORS[1], "--"),
    ("c_mockEE", r"analytic, mock mean $C_\ell^{EE}$", CATEGORICAL_COLORS[2], ":"),
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mixing-test", required=True)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    with open(a.mixing_test) as fh:
        d = json.load(fh)

    ell = np.array(d["ell_eff"])
    emp = np.array(d["sigma_BB_empirical"])
    cut_lo, cut_hi = d["reporting_cut"]

    plt.style.use(PAPER_MPLSTYLE)
    fig, (ax, rx) = plt.subplots(
        2, 1, figsize=(FIG_WIDTH_FULL, 5.0), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.08},
    )

    for axis in (ax, rx):
        axis.axvspan(cut_lo, cut_hi, color="0.85", alpha=0.45, lw=0, zorder=0)

    ax.plot(ell, emp, color="0.15", lw=1.6, marker="o", ms=3.2, zorder=5,
            label=rf"empirical, {d['n_realisations']} mocks")
    for key, label, colour, ls in VARIANTS:
        ax.plot(ell, np.array(d["variants"][key]["sigma_BB"]), color=colour,
                ls=ls, lw=1.4, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"$\sigma(C_\ell^{BB})$")
    ax.legend(frameon=False, fontsize="small", loc="lower left")
    ax.set_title(
        r"E-to-B mixing dominates the low-$\ell$ analytic $C_\ell^{BB}$ covariance",
        fontsize="medium",
    )

    rx.axhline(1.0, color="0.15", lw=0.9, zorder=4)
    for key, label, colour, ls in VARIANTS:
        rx.plot(ell, np.array(d["variants"][key]["ratio_emp_over_ana"]),
                color=colour, ls=ls, lw=1.4)
    rx.set_xscale("log")
    rx.set_xlim(8, 2048)
    rx.set_ylim(0.0, 3.6)
    rx.set_xlabel(r"$\ell$")
    rx.set_ylabel(r"$\sigma_{\rm emp} / \sigma_{\rm ana}$")

    fig.savefig(a.out_file, bbox_inches="tight", dpi=200)
    print(f"wrote {a.out_file}")


if __name__ == "__main__":
    main()
