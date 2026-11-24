"""Combined GLASS-mock validation figure for the B-mode paper.

Two rows over the three B-mode statistics.  --range fiducial (the paper
figure) shows each statistic on its own scale cuts; --range full shows every
reporting bin, with the excluded scales shaded and the COSEBIs both at the
fiducial cuts (filled) and over all scales (open).  Row 1 is the across-realization mean in units of the analytic
(data) sigma, with
sigma_emp/(sigma_ana sqrt(N)) error bars, so the y-axis reads as "fraction
of a data error bar".  Row 2 checks that the analytic covariance reproduces
the mock scatter, as sigma_ana / sigma_emp so that a value below one reads "the analytic covariance is tighter than the mock
scatter", with the +-1/sqrt(2N) sampling uncertainty of that ratio on each point.

Columns: COSEBI B_n | xi_+^B and xi_-^B | C_l^BB.  Colour encodes the
estimator family, as in the data-vector figures: purple for COSEBIs (light
tint for all scales), crimson for xi_+^B and orange for xi_-^B, teal for
C_l^BB; the second series of a panel is also the open marker.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import sys  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# Resolved from this file, as the other workflow scripts do, so the figure
# builds from any working directory.
SRC = Path(__file__).resolve().parents[2]
PAPER_MPLSTYLE = str(SRC / "workflow" / "matplotlib_config" / "paper.mplstyle")

# The "squareroot" x scale is registered as a side effect of importing the
# paper's plotting helpers; the harmonic panel uses it so the power-0.5
# bandpowers land evenly, as in the C_l data-vector figure.
sys.path.insert(0, str(SRC / "papers" / "bmodes" / "scripts"))
from plotting_utils import (  # noqa: E402
    COLOR_CL_BB,
    COLOR_COSEBIS,
    COLOR_COSEBIS_ALL,
    COLOR_XI_B,
    COLOR_XI_B_MINUS,
    DATA_VECTOR_MARK,
    EXCLUDED_SHADE,
    FIG_WIDTH_FULL,
)

# Fiducial windows: config-space cut in arcmin, harmonic cut in multipole.
THETA_CUT = (12.0, 83.0)
ELL_CUT = (300.0, 1600.0)

ROW1_YLIM = (-0.25, 0.25)
ROW2_YLIM = {"fiducial": (0.8, 1.2), "full": (0.7, 1.5)}
# The lowest C_l band's sigma ratio (~1.8) would set the whole ratio row's
# range on its own; it is dropped from that panel (kept in the mean panel).
RATIO_SKIP_LOWEST_CL = 1
ELL_TICKS = [100, 400, 900, 1600]

# Drawn at the A&A text width, the width the paper prints it at, so the
# style's font sizes are the ones the reader gets.
FIG_HEIGHT = 2.95

# Ranges are named, never spelled out in arcminutes or multipoles: the
# words carry across both spaces and the numbers are defined in the text.

MARK = DATA_VECTOR_MARK


def load_full_cl(cl_dir: Path, mock_ids) -> np.ndarray:
    """Stack the per-mock 32-band C_l^BB vectors.

    The campaign npz stores only the 17 in-cut bands, so the full range is
    read from the same per-mock npy files mock_campaign_analysis.py reads
    (GLASS convention: row 0 is ell_eff, row 4 is C_l^BB).
    """
    vectors = []
    for mock_id in mock_ids:
        path = cl_dir / f"cl_glass_mock_{mock_id}_4096.npy"
        array = np.load(path, allow_pickle=False)
        vectors.append(array[4])
    return np.asarray(vectors)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-evidence", required=True)
    parser.add_argument("--per-mock-statistics", required=True)
    parser.add_argument(
        "--cl-dir",
        default=None,
        help="Per-mock external C_l npy directory, needed only by --range "
        "full. Defaults to the directory the campaign recorded in its "
        "provenance, which is the one the stack was actually built from.",
    )
    parser.add_argument(
        "--range",
        dest="scale_range",
        choices=("fiducial", "full"),
        default="fiducial",
        help="Fiducial cuts only (the paper figure), or every reporting bin.",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    plt.style.use(PAPER_MPLSTYLE)

    data = np.load(args.per_mock_statistics, allow_pickle=False)
    provenance = json.loads(str(data["provenance"]))
    cl_dir = Path(args.cl_dir or provenance["ensembles"]["cl_dir"])
    with open(args.campaign_evidence) as stream:
        evidence = json.load(stream)["sigma_ratios"]

    n_real = data["Bn"].shape[0]
    ratio_err = 1.0 / np.sqrt(2.0 * n_real)
    theta = np.asarray(data["theta"])
    ell = np.asarray(data["ell"])
    modes = np.arange(1, data["Bn"].shape[1] + 1)
    cl_full = None
    if args.scale_range == "full":
        cl_full = load_full_cl(cl_dir, [str(v) for v in data["mock_ids"]])
    # The stored in-cut stack is the same bytes, masked: the guard catches a
    # wrong cl_dir, which otherwise only shows up as a distorted panel.
    if cl_full is not None and not np.allclose(
        cl_full[:, np.asarray(data["cl_mask"])], data["cl_BB"]
    ):
        raise SystemExit(
            f"{cl_dir} does not reproduce the campaign's in-cut C_l stack; "
            "the per-mock C_l directory is wrong."
        )

    def summarize(vectors, sigma):
        """Mean in units of sigma_ana, its error, and sigma_ana / sigma_emp."""
        vectors = np.asarray(vectors)
        sigma = np.asarray(sigma)
        scatter = vectors.std(axis=0, ddof=1)
        return (
            vectors.mean(axis=0) / sigma,
            scatter / sigma / np.sqrt(n_real),
            sigma / scatter,
        )

    fiducial = args.scale_range == "fiducial"
    # column, x, (mean, err, ratio), marker, filled, colour, label.  Within a
    # panel the fiducial series (xi_+^B, COSEBIs at the fiducial cuts) is the
    # filled circle in the family hue; the other one (xi_-^B, COSEBIs over
    # all scales) the open diamond in its own hue.
    if fiducial:
        theta_cut = np.asarray(data["theta_cut"])
        series = [
            (0, modes, summarize(data["Bn"], data["sigma_cosebis"]),
             "o", True, COLOR_COSEBIS, None),
            (1, theta_cut * 0.975,
             summarize(data["xip_B_cut"], data["sigma_xip_cut"]),
             "o", True, COLOR_XI_B, r"$\xi_+^B$"),
            (1, theta_cut * 1.025,
             summarize(data["xim_B_cut"], data["sigma_xim_cut"]),
             "D", False, COLOR_XI_B_MINUS, r"$\xi_-^B$"),
            (2, np.asarray(data["ell_cut"]),
             summarize(data["cl_BB"], data["sigma_cl_cut"]),
             "o", True, COLOR_CL_BB, None),
        ]
    else:
        series = [
            (0, modes - 0.18, summarize(data["Bn"], data["sigma_cosebis"]),
             "o", True, COLOR_COSEBIS, "fiducial cuts"),
            (0, modes + 0.18,
             summarize(data["Bn_full"], data["sigma_cosebis_full"]),
             "D", False, COLOR_COSEBIS_ALL, "all scales"),
            (1, theta * 0.965,
             summarize(data["xip_B_full"], data["sigma_xip_full"]),
             "o", True, COLOR_XI_B, r"$\xi_+^B$"),
            (1, theta * 1.035,
             summarize(data["xim_B_full"], data["sigma_xim_full"]),
             "D", False, COLOR_XI_B_MINUS, r"$\xi_-^B$"),
            (2, ell, summarize(cl_full, data["sigma_cl_full"]),
             "o", True, COLOR_CL_BB, None),
        ]

    fig, axes = plt.subplots(
        2, 3, figsize=(FIG_WIDTH_FULL, FIG_HEIGHT), sharey="row", sharex="col",
    )

    reported = {}
    for column, x, (mean, err, ratio), marker, filled, color, label in series:
        face = color if filled else "none"
        style = dict(MARK, color=color, ecolor=color, marker=marker,
                     markerfacecolor=face, markeredgecolor=color)
        axes[0][column].errorbar(x, mean, yerr=err, label=label, **style)
        skip = RATIO_SKIP_LOWEST_CL if (column == 2 and not fiducial) else 0
        axes[1][column].errorbar(x[skip:], ratio[skip:], yerr=ratio_err, **style)
        reported[label or f"column {column}"] = (float(np.max(np.abs(mean))),
                           float(ratio.min()), float(ratio.max()))

    # Excluded scales shaded, as in the paper's data-vector figures.  The
    # fiducial figure has nothing outside its cuts, so it carries no shading.
    axes[1][1].set_xscale("log")
    axes[1][2].set_xscale("squareroot")
    if fiducial:
        axes[1][1].set_xlim(11.0, 95.0)
        axes[1][2].set_xlim(ELL_CUT[0] * 0.9, ELL_CUT[1] * 1.05)
    else:
        axes[1][1].set_xlim(0.9, 290.0)
        axes[1][2].set_xlim(8.0, 2100.0)
        for row in (0, 1):
            for column, cut in ((1, THETA_CUT), (2, ELL_CUT)):
                ax = axes[row][column]
                low, high = ax.get_xlim()
                ax.axvspan(low, cut[0], **EXCLUDED_SHADE)
                ax.axvspan(cut[1], high, **EXCLUDED_SHADE)

    for ax in axes[0]:
        ax.axhline(0.0, color="black", ls="--", lw=0.8, alpha=0.7, zorder=1)
        ax.set_ylim(ROW1_YLIM)
    for ax in axes[1]:
        ax.axhline(1.0, color="black", ls="--", lw=0.8, alpha=0.7, zorder=1)
        ax.set_ylim(ROW2_YLIM[args.scale_range])

    for ax, text in zip(axes[0], ("COSEBIs", r"$\xi^B_\pm$", r"$C_\ell^{BB}$")):
        ax.text(
            0.5, 0.95, text, transform=ax.transAxes, ha="center", va="top",
            fontsize=plt.rcParams["axes.labelsize"],
        )

    axes[1][0].set_xlabel(r"COSEBI mode $n$")
    axes[1][1].set_xlabel(r"$\theta$ [arcmin]")
    axes[1][2].set_xlabel(r"$\ell$")
    axes[0][0].set_ylabel(r"$\langle B \rangle / \sigma^{\rm ana}$")
    axes[1][0].set_ylabel(r"$\sigma^{\rm ana} / \sigma^{\rm emp}$")
    axes[1][0].set_xticks([1, 5, 10, 15, 20])
    theta_ticks = [20, 50, 80] if fiducial else [1, 10, 100]
    axes[1][1].set_xticks(theta_ticks)
    axes[1][1].set_xticklabels([str(value) for value in theta_ticks])
    axes[1][1].xaxis.set_minor_formatter(ticker.NullFormatter())
    # Only ticks the panel actually spans: a label pinned to the frame at a
    # value outside the limits misreads as the axis starting there.
    low, high = axes[1][2].get_xlim()
    axes[1][2].set_xticks([t for t in ELL_TICKS if low <= t <= high])

    def handle(marker, filled, color, label):
        return Line2D([], [], ls="none", marker=marker, color=color,
                      markerfacecolor=color if filled else "none",
                      markeredgecolor=color,
                      markeredgewidth=MARK["markeredgewidth"],
                      markersize=MARK["markersize"], label=label)

    # The statistic label owns the top-left of the top pure-mode panel and
    # the data crowd its top-right; the horizontal xi_pm legend takes the
    # top-right of the ratio panel below, which is empty in both views.
    legend_kw = {"frameon": False, "handletextpad": 0.3, "borderpad": 0.2,
                 "columnspacing": 0.8, "borderaxespad": 0.5, "ncol": 2}
    axes[1][1].legend(
        handles=[handle("o", True, COLOR_XI_B, r"$\xi_+^B$"),
                 handle("D", False, COLOR_XI_B_MINUS, r"$\xi_-^B$")],
        loc="upper right", **legend_kw)
    if not fiducial:
        # Same placement for the COSEBI pair: the ratio panel's top-right is
        # empty (the low-n modes sit below one) and the top panel's is not.
        axes[1][0].legend(
            handles=[handle("o", True, COLOR_COSEBIS, "fiducial cuts"),
                     handle("D", False, COLOR_COSEBIS_ALL, "all scales")],
            loc="upper right", **legend_kw)

    fig.align_ylabels(axes[:, 0])
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    print(f"Saved {out} ({FIG_WIDTH_FULL} x {FIG_HEIGHT} in)")
    for label, (max_mean, low, high) in reported.items():
        print(f"{label}: max |mean|/sigma_ana = {max_mean:.3f}; "
              f"sigma_ana/sigma_emp {low:.3f}-{high:.3f}")
    if cl_full is not None:
        print(f"full-range C_l read from {cl_dir}")
    # The evidence file is the reference for the ranges it does cover.
    for key, series_key, sigma_key in (
        ("COSEBI_B_modes_1_to_20", "Bn", "sigma_cosebis"),
        ("C_ell_BB_17_cut_bands", "cl_BB", "sigma_cl_cut"),
    ):
        ours = np.asarray(data[series_key]).std(axis=0, ddof=1) / data[sigma_key]
        theirs = np.asarray(evidence[key]["ratio"])
        print(f"cross-check {key}: max |delta| = "
              f"{np.max(np.abs(ours - theirs)):.2e}")
    plt.close(fig)


if __name__ == "__main__":
    main()
