"""The scan's other half: how NOISY each variant is, not just how biased.

``pure_eb_decision_comparison_figure`` reads every panel in units of the
variant's own analytic sigma -- mean/sigma_ana, sigma_emp/sigma_ana -- which is
exactly the right frame for "is this variant biased, is its covariance
calibrated" and exactly the wrong one for "which variant would you rather
measure with".  A ratio of two numbers that both grow says nothing about
either.  A variant can sit at mean/sigma = 0.1 with an error bar twice
anyone else's and read as the best row in that figure.

So this one plots the scatter itself.  Per variant, per reporting bin, the
across-realisation standard deviation of the mock pure-B vector -- the mocks
carry no B power, so that sigma IS the method's noise -- multiplied by theta
in arcmin, which flattens the curves enough that seven of them can share one
pair of axes.  Units are therefore arcmin (theta) times the dimensionless
xi^B, quoted in 1e-6 arcmin.

Two panels, xi_+^B and xi_-^B, sharing the theta axis.  Seven lines, one per
variant, in the data-viz reference categorical order with a distinct dash
pattern each, the adopted variant heavy and solid and the published one
dashed.  The adopted variant's *analytic* theta.sigma -- the propagated
covariance's diagonal, from its covariance bundle -- is drawn thin and grey
behind them, so the empirical spread of every variant can be read against the
one error bar the paper actually quotes.  The adopted [12, 83] arcmin window
is shaded.

    python pure_eb_scan_scatter_figure.py --mock-modes M.npz \
        --covariance-bundle B.tar.gz --out-file F.png
"""

from __future__ import annotations

import argparse
import tarfile
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mock_campaign_analysis import (
    PURE_THETA_BY_CONVENTION,
    _load_pure_covariance,
)
from plotting_utils import FIG_WIDTH_FULL, PAPER_MPLSTYLE
from pure_eb_scan_variants import ADOPTED, BY_NAME, PUBLISHED, add_variants_argument, resolve

plt.style.use(PAPER_MPLSTYLE)

#: Slots 1-7 of the data-viz reference categorical palette, in its documented
#: order (blue, orange, aqua, yellow, magenta, green, violet).  Seven line
#: series is past the three-slot all-pairs cap, so the adjacent-pair gate is
#: the one in force -- which this order clears by construction -- and every
#: series carries a dash pattern as secondary encoding on top of hue.
SERIES_COLORS = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
)
#: Dash patterns, one per series slot.  Slot 1 is the adopted variant's solid
#: line and slot 2 the published variant's dash; the rest are distinguishable
#: at 2pt on a 3.5in panel.
SERIES_DASHES = (
    "solid",
    (0, (5, 2)),
    (0, (1, 1.4)),
    (0, (4, 1.4, 1, 1.4)),
    (0, (3, 1.2, 1, 1.2, 1, 1.2)),
    (0, (7, 2, 2, 2)),
    (0, (2, 1.6)),
)
#: Adopted first, published second, the rest in scan-table order.  Colour
#: follows the variant, not its rank, so this mapping is fixed regardless of
#: which variants a run selects.
SERIES_ORDER = (ADOPTED, PUBLISHED)
SHADE = dict(color="0.90", zorder=0)
#: theta.sigma runs at a few times 1e-6; quoting the axis in those units keeps
#: the tick labels to one significant figure instead of an exponent per tick.
Y_SCALE = 1e-6
PANELS = ((r"$\xi_+^B$", "xip_B", "sigma_xip_full"), (r"$\xi_-^B$", "xim_B", "sigma_xim_full"))


def _ordered(variants):
    """Scan variants in drawing order: adopted, published, then table order."""

    lead = [BY_NAME[name] for name in SERIES_ORDER if name in {v.name for v in variants}]
    rest = [variant for variant in variants if variant.name not in SERIES_ORDER]
    return lead + rest


def _theta_sigma(modes, variant, key):
    """theta [arcmin] times the across-mock standard deviation of a mode block."""

    theta = np.asarray(modes[f"{variant.name}__theta"])[0]
    values = np.asarray(modes[f"{variant.name}__{key}"])
    return theta, theta * values.std(axis=0, ddof=1)


def _analytic(bundle_path, reference, theta):
    """theta.sigma_ana for both mode blocks, from one covariance bundle."""

    with tempfile.TemporaryDirectory() as tmp:
        staging = Path(tmp)
        with tarfile.open(bundle_path, "r:gz") as bundle:
            bundle.extractall(staging)
        members = sorted(staging.glob(f"*__{reference}.npz"))
        if len(members) != 1:
            raise ValueError(
                f"expected exactly one *__{reference}.npz in {bundle_path}, "
                f"found {[path.name for path in members]}"
            )
        covariance = _load_pure_covariance(members[0], theta)
    return {key: theta * covariance[key] for _, _, key in PANELS}


def _table(rows, adopted_name):
    """Window-median theta.sigma_emp per variant, relative to the adopted one."""

    reference = rows[adopted_name]
    width = max(len(name) for name in rows)
    lines = [
        f"{'variant'.ljust(width)}  "
        f"{'xi+^B':>10} {'/adopted':>9}  {'xi-^B':>10} {'/adopted':>9}"
    ]
    for name, values in rows.items():
        cells = []
        for key in ("xip_B", "xim_B"):
            cells.append(f"{values[key] / Y_SCALE:>10.3f}")
            cells.append(f"{values[key] / reference[key]:>9.3f}")
        lines.append(f"{name.ljust(width)}  {cells[0]} {cells[1]}  {cells[2]} {cells[3]}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mock-modes", required=True)
    parser.add_argument(
        "--covariance-bundle",
        required=True,
        help=(
            "The ADOPTED variant's covariance bundle from "
            "pure_eb_scan_covariances; its diagonal is the grey analytic "
            "reference the empirical curves are read against."
        ),
    )
    parser.add_argument(
        "--reference-covariance",
        choices=("data", "mockmatched"),
        default="mockmatched",
        help=(
            "Which of the bundle's two covariances the analytic reference "
            "comes from. Default 'mockmatched': the curves it sits behind are "
            "the scatter of the MOCKS, so the like-for-like analytic sigma is "
            "the same propagation at the mocks' noise level."
        ),
    )
    parser.add_argument(
        "--scale-cut-convention",
        choices=tuple(PURE_THETA_BY_CONVENTION),
        default="grid_bin_edges",
    )
    parser.add_argument("--out-file", required=True)
    add_variants_argument(parser)
    args = parser.parse_args()

    variants = _ordered(resolve(args.variants))
    window = PURE_THETA_BY_CONVENTION[args.scale_cut_convention]

    figure, axes = plt.subplots(
        1, 2, figsize=(FIG_WIDTH_FULL, 0.42 * FIG_WIDTH_FULL), sharex=True
    )
    medians: dict[str, dict[str, float]] = {}
    with np.load(args.mock_modes, allow_pickle=False) as modes:
        available = {str(value) for value in np.asarray(modes["variants"])}
        missing = [variant.name for variant in variants if variant.name not in available]
        if missing:
            raise ValueError(
                f"{missing} absent from the mock pure-mode scan "
                f"({sorted(available)}); the two stages must run the same "
                "variant list."
            )
        for index, variant in enumerate(variants):
            color = SERIES_COLORS[index % len(SERIES_COLORS)]
            dashes = SERIES_DASHES[index % len(SERIES_DASHES)]
            medians[variant.name] = {}
            for axis, (_, key, _) in zip(axes, PANELS):
                theta, scatter = _theta_sigma(modes, variant, key)
                inside = (theta >= window[0]) & (theta <= window[1])
                medians[variant.name][key] = float(np.median(scatter[inside]))
                axis.plot(
                    theta,
                    scatter / Y_SCALE,
                    color=color,
                    ls=dashes,
                    lw=1.6 if variant.name == ADOPTED else 1.0,
                    label=variant.label,
                    zorder=3 if variant.name == ADOPTED else 2,
                )
        reference_theta = np.asarray(modes[f"{ADOPTED}__theta"])[0]

    analytic = _analytic(args.covariance_bundle, args.reference_covariance, reference_theta)
    for axis, (title, _, sigma_key) in zip(axes, PANELS):
        axis.plot(
            reference_theta,
            analytic[sigma_key] / Y_SCALE,
            color="0.45",
            lw=0.9,
            zorder=1,
            label=r"adopted $\sigma_{\rm ana}$",
        )
        axis.axvspan(*window, **SHADE)
        axis.set_xscale("log")
        # Log in y as well: the spread between variants is a factor of 6.5,
        # and on a linear axis the four variants that agree to a percent
        # collapse onto one another at the bottom of the frame.
        axis.set_yscale("log")
        axis.set_yticks([4, 5, 6, 8, 10, 15, 20, 30])
        axis.set_yticklabels(["4", "5", "6", "8", "10", "15", "20", "30"])
        axis.minorticks_off()
        axis.set_xlabel(r"$\theta$ [arcmin]")
        axis.set_title(title)
    axes[0].set_ylabel(r"$\theta\,\sigma(\theta)$ [$10^{-6}$ arcmin]")

    # The legend goes below both panels rather than inside either: at seven
    # entries plus the reference it is taller than a panel, and the curves
    # leave no corner empty once the analytic line is in.
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=3,
        frameon=False,
        fontsize=7,
        handlelength=3.2,
    )
    figure.tight_layout()

    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out_path, bbox_inches="tight")
    print(_table(medians, ADOPTED))


if __name__ == "__main__":
    main()
