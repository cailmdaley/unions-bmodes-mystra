"""The pure-E/B decision scan, as one paper figure.

One row per decision variant, over two shared columns:

* left   -- the across-realisation mean xi_+^B / xi_-^B in units of the
            analytic (data) sigma, with sigma_emp/(sigma_ana sqrt(N)) error
            bars.  The mocks carry no B power, so any excursion from zero is
            that variant's own rectification bias, and the y-axis reads as
            "fraction of a data error bar";
* right  -- the empirical-over-analytic sigma ratio, with the 1/sqrt(2(N-1))
            sampling error of an N-realisation standard deviation.  Says
            whether the variant's propagated covariance has the right
            amplitude.

Variants read down the page, which leaves room for the labels and the numbers.
Each row carries its identity as a left-margin label -- the
adopted variant bold, the published one italic -- and its numbers in a third,
axis-free text column at the right, aligned as a small table so the mock
verdict, the KS p-value and the data PTE can be read down the page as three
columns rather than hunted for inside the panels.

Both panels are scatter with error bars rather than connected lines -- these
are independent per-bin estimates -- and the two series are dodged in theta so
their markers do not overlap.  The adopted scale cuts are shaded.

Everything is a function of the scan statistics json alone, so re-cutting the
figure never re-runs a transform.  A&A double-column width, the paper
matplotlib style, no suptitle.

    python pure_eb_decision_comparison_figure.py --stats S.json --out-file F.png
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plotting_utils import (
    CATEGORICAL_COLORS,
    FIG_WIDTH_FULL,
    PAPER_MPLSTYLE,
    SCATTER_ERRORBAR,
)

plt.style.use(PAPER_MPLSTYLE)
# The margins here are chosen, not derived: the left one holds the row labels
# and the right one the numbers table, neither of which tight_layout can see.
plt.rcParams["figure.autolayout"] = False

SHADE = dict(color="0.90", zorder=0)
STATS_FONTSIZE = 7
LABEL_FONTSIZE = 6.5
# Marker shape repeats the colour distinction for greyscale and CVD readers.
SERIES = ((r"$\xi_+^B$", "xi+^B", "o"), (r"$\xi_-^B$", "xi-^B", "s"))
# The two series are estimated in the same theta bins; without a nudge their
# markers and error bars land on top of each other.  +/-2.5%, far smaller than
# a bin, multiplicative because the axis is logarithmic.
DODGE = 0.025
# Columns of the numbers table, in axes coordinates of its (invisible) axes:
# the row name is left-aligned at 0, the three values right-aligned under
# their headers.
STAT_LABEL_X = -0.23
STAT_VALUE_X = (0.52, 0.76, 1.00)
STAT_ROW_Y = (0.70, 0.50, 0.30)
# Grid columns the three panel columns occupy; 1 and 3 are bare spacers.
COLUMN_SLOTS = (0, 2, 4)


def _p(value):
    """A p-value in the fewest characters that still says what it is."""

    value = float(value)
    return f"{value:.0e}".replace("e-0", "e-") if value < 1e-3 else f"{value:.2g}"


def _panels(axes, entry, theta_cut, colors):
    """The two data panels of one variant's row."""

    mean_axis, ratio_axis = axes
    for index, ((label, name, marker), color) in enumerate(zip(SERIES, colors)):
        profile = entry["profiles"][name]
        theta = np.asarray(profile["theta"], dtype=float)
        theta = theta * (1.0 + DODGE * (2.0 * index / (len(SERIES) - 1) - 1.0))
        marks = dict(
            marker=marker, color=color, ecolor=color, label=label, **SCATTER_ERRORBAR
        )
        mean_axis.errorbar(
            theta,
            np.asarray(profile["mean_over_sigma_ana"], dtype=float),
            yerr=np.asarray(profile["stderr_over_sigma_ana"], dtype=float),
            **marks,
        )
        ratio_axis.errorbar(
            theta,
            np.asarray(profile["sigma_ratio"], dtype=float),
            yerr=profile["sigma_ratio_error"],
            **marks,
        )
    for axis, reference in ((mean_axis, 0.0), (ratio_axis, 1.0)):
        axis.axhline(reference, color="black", ls="--", lw=0.7, zorder=1)
        axis.axvspan(*theta_cut, **SHADE)
        axis.set_xscale("log")


def _numbers(axis, entry, header):
    """One variant's numbers, as a three-row table in an axis-free column."""

    axis.set_axis_off()
    per_mock = entry["per_mock"]
    data = entry["data_pte"]
    rows = (
        (
            r"mock $\langle B\rangle$ PTE",
            (
                _p(per_mock["xi+^B"]["mean_vector"]["pte"]),
                _p(per_mock["xi-^B"]["mean_vector"]["pte"]),
                "",
            ),
        ),
        (
            "KS",
            (_p(per_mock["xi+^B"]["ks_pvalue"]), _p(per_mock["xi-^B"]["ks_pvalue"]), ""),
        ),
        (
            "data PTE",
            (_p(data["xi+^B"]), _p(data["xi-^B"]), _p(data["xi_tot^B"])),
        ),
    )
    if header:
        for x, title in zip(STAT_VALUE_X, (r"$\xi_+$", r"$\xi_-$", "tot")):
            axis.text(
                x, 0.90, title, transform=axis.transAxes, fontsize=STATS_FONTSIZE,
                ha="right", va="center", color="0.35",
            )
    for y, (name, values) in zip(STAT_ROW_Y, rows):
        axis.text(
            STAT_LABEL_X, y, name, transform=axis.transAxes,
            fontsize=STATS_FONTSIZE, ha="left", va="center",
        )
        for x, value in zip(STAT_VALUE_X, values):
            axis.text(
                x, y, value, transform=axis.transAxes, fontsize=STATS_FONTSIZE,
                ha="right", va="center", family="monospace",
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args()

    payload = json.loads(Path(args.stats).read_text())
    order = payload["variant_order"]
    theta_cut = tuple(payload["theta_cut"])
    colors = CATEGORICAL_COLORS[: len(SERIES)]

    n = len(order)
    fig = plt.figure(figsize=(FIG_WIDTH_FULL, 1.55 * n + 0.9))
    # Explicit spacer columns rather than a uniform wspace: the gap before the
    # ratio panel has to clear its y tick labels, the gap before the numbers
    # only has to separate two blocks, and they are not the same width.
    grid = fig.add_gridspec(
        n, 5, width_ratios=[1.0, 0.32, 1.0, 0.30, 1.20], wspace=0.0, hspace=0.14,
    )
    axes = np.empty((n, 3), dtype=object)
    for row in range(n):
        for column, slot in enumerate(COLUMN_SLOTS):
            # One y-range per column, so a bias that is large in one variant
            # reads as large next to the others rather than being rescaled
            # away; x is shared everywhere.
            axes[row, column] = fig.add_subplot(
                grid[row, slot],
                sharex=axes[0, column] if row else None,
                sharey=axes[0, column] if row else None,
            )
            if row < n - 1:
                axes[row, column].tick_params(labelbottom=False)

    for index, name in enumerate(order):
        entry = payload["variants"][name]
        _panels(axes[index, :2], entry, theta_cut, colors)
        _numbers(axes[index, 2], entry, header=index == 0)
        # The row label carries the tag in its own weight, so naming the
        # adopted and published variants costs no extra line of text.
        # Broken onto its own lines: the margin is a column, not a caption.
        lines = re.split(r",\s*|\s+(?=\()", entry["label"])
        weight, style = "normal", "normal"
        if name == payload["adopted"]:
            lines.append("[adopted]")
            weight = "bold"
        elif name == payload["published"]:
            lines.append("[published]")
            style = "italic"
        axes[index, 0].text(
            -0.72, 0.5, "\n".join(lines), transform=axes[index, 0].transAxes,
            fontsize=LABEL_FONTSIZE, ha="center", va="center",
            fontweight=weight, fontstyle=style,
        )

    axes[0, 0].set_title(r"$\langle\xi^B\rangle/\sigma^\mathrm{ana}$", fontsize=9)
    axes[0, 1].set_title(r"$\sigma^\mathrm{emp}/\sigma^\mathrm{ana}$", fontsize=9)
    for axis in axes[-1, :2]:
        axis.set_xlabel(r"$\theta$ [arcmin]")
    axes[0, 0].legend(fontsize=6.5, loc="lower left", ncol=2, handletextpad=0.3,
                      columnspacing=1.0, borderpad=0.3, framealpha=0.9)

    fig.subplots_adjust(left=0.320, right=0.985, top=0.952, bottom=0.062)
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
