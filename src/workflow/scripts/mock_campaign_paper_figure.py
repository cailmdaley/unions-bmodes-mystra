"""The GLASS mock campaign, as one paper figure.

Two rows over the same three columns (COSEBI modes, pure ξ^B, C_ℓ^BB):

* top    — the across-realisation mean B vector in units of the analytic
           (data) σ, with σ_emp/(σ_ana √N) error bars.  Reads as "fraction of
           a data error bar";
* bottom — the empirical-over-analytic σ ratio, with the 1/sqrt(2(N-1))
           sampling error of an N-realisation standard deviation.

Both rows are functions of the campaign's per-realisation statistics npz
alone, and both reuse ``mock_campaign_analysis``'s loader and cut constants
rather than restating them: the two figures this composes
(``mock_mean_vector_figure`` and ``mock_sigma_ratio_figure``) stay the
working versions, and this is their paper-format twin — no suptitle, A&A
double-column width, the paper matplotlib style.

The adopted scale cuts are shaded in the θ column.  The C_ℓ^BB column carries
no band: every bandpower the campaign stores for it already lies inside
300 < ℓ < 1600, so a band would grey the whole panel.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plotting_utils import (
    CATEGORICAL_COLORS,
    FIG_WIDTH_FULL,
    PAPER_MPLSTYLE,
    SCATTER_ERRORBAR,
)
from mock_campaign_analysis import (
    PURE_THETA_MAX,
    PURE_THETA_MIN,
    read_statistics,
)

plt.style.use(PAPER_MPLSTYLE)

SHADE = dict(color="0.90", zorder=0)
ANNOTATION_FONTSIZE = 8
ANNOTATION_BOX = dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5)


def _chi2_pte_of_mean(vectors: np.ndarray) -> tuple[float, float, int]:
    """χ² of the mean against the empirical covariance of the mean."""

    from scipy import stats

    n_real, n_bins = vectors.shape
    mean = vectors.mean(axis=0)
    cov_mean = np.cov(vectors, rowvar=False, ddof=1) / n_real
    hartlap = (n_real - n_bins - 2) / (n_real - 1)
    chi2 = float(mean @ np.linalg.inv(cov_mean) @ mean * hartlap)
    return chi2, float(stats.chi2.sf(chi2, n_bins)), n_bins


def _dodge(x, index, count, log):
    """Nudge paired series apart at a shared abscissa.

    ξ_+^B and ξ_-^B are estimated in the *same* θ bins, so their markers and
    error bars would sit on top of each other.  The offset is ±2.5% — visibly
    smaller than a bin — and applied multiplicatively on a log axis so it is
    the same nudge everywhere.
    """

    if count < 2:
        return x
    shift = 0.025 * (2.0 * index / (count - 1) - 1.0)
    return x * (1.0 + shift) if log else x + shift * np.ptp(x) / len(x)


def _mean_panel(ax, x, vectors, sigma_ana, label, marker, color):
    n_real = vectors.shape[0]
    mean = vectors.mean(axis=0) / sigma_ana
    err = vectors.std(axis=0, ddof=1) / sigma_ana / np.sqrt(n_real)
    ax.errorbar(
        x,
        mean,
        yerr=err,
        marker=marker,
        color=color,
        ecolor=color,
        label=label,
        **SCATTER_ERRORBAR,
    )
    chi2, pte, dof = _chi2_pte_of_mean(vectors)
    return rf"{label}: $\chi^2$={chi2:.1f}/{dof}, $p$={pte:.2g}"


def _ratio_panel(ax, x, vectors, sigma_ana, label, marker, color, rel_err):
    ratio = vectors.std(axis=0, ddof=1) / sigma_ana
    ax.errorbar(
        x,
        ratio,
        yerr=ratio * rel_err,
        marker=marker,
        color=color,
        ecolor=color,
        label=label,
        **SCATTER_ERRORBAR,
    )
    return ratio


def _symmetric_limits(values, centre):
    finite = np.asarray(values)[np.isfinite(values)]
    deviation = np.max(np.abs(finite - centre))
    return (centre - 1.05 * deviation, centre + 1.05 * deviation)


def build_figure(payload, out_path: Path) -> None:
    n_real = payload["Bn"].shape[0]
    rel_err = 1.0 / np.sqrt(2.0 * (n_real - 1))
    modes = np.arange(1, 21)
    theta = np.asarray(payload["theta"])
    ell = np.asarray(payload["ell_cut"])

    # Fixed categorical order, never cycled: slot 1 carries the single-series
    # panels and ξ_+^B, slot 2 the only series that ever shares a panel with
    # it.  Marker shape repeats the distinction for greyscale and CVD readers.
    palette = CATEGORICAL_COLORS

    columns = (
        (
            [("$B_n$", modes, payload["Bn"], payload["sigma_cosebis"], "o", palette[0])],
            "COSEBI mode $n$",
            None,
            None,
        ),
        (
            [
                (
                    r"$\xi_+^B$",
                    theta,
                    payload["xip_B_full"],
                    payload["sigma_xip_full"],
                    "o",
                    palette[0],
                ),
                (
                    r"$\xi_-^B$",
                    theta,
                    payload["xim_B_full"],
                    payload["sigma_xim_full"],
                    "s",
                    palette[1],
                ),
            ],
            r"$\theta$ [arcmin]",
            "log",
            (PURE_THETA_MIN, PURE_THETA_MAX),
        ),
        (
            [
                (
                    r"$C_\ell^{BB}$",
                    ell,
                    payload["cl_BB"],
                    payload["sigma_cl_cut"],
                    "o",
                    palette[0],
                )
            ],
            r"$\ell$",
            "squareroot",
            None,
        ),
    )

    fig, axes = plt.subplots(
        2, 3, figsize=(FIG_WIDTH_FULL, 4.4), sharex="col", squeeze=False
    )
    ratio_content = []
    mean_content = []
    for column, (series, xlabel, xscale, shaded) in enumerate(columns):
        top, bottom = axes[0, column], axes[1, column]
        annotations = []
        for index, (label, x, vectors, sigma_ana, marker, colour) in enumerate(series):
            x = _dodge(np.asarray(x, dtype=float), index, len(series), xscale == "log")
            vectors = np.asarray(vectors)
            sigma_ana = np.asarray(sigma_ana)
            annotations.append(
                _mean_panel(top, x, vectors, sigma_ana, label, marker, colour)
            )
            ratio = _ratio_panel(
                bottom, x, vectors, sigma_ana, label, marker, colour, rel_err
            )
            ratio_content.append(ratio * (1.0 + rel_err))
            ratio_content.append(ratio * (1.0 - rel_err))
            mean = vectors.mean(axis=0) / sigma_ana
            spread = vectors.std(axis=0, ddof=1) / sigma_ana / np.sqrt(n_real)
            mean_content.extend([mean + spread, mean - spread])
        top.axhline(0.0, color="black", ls="--", lw=0.8)
        bottom.axhline(1.0, color="black", ls="--", lw=0.8)
        bottom.set_xlabel(xlabel)
        for ax in (top, bottom):
            if xscale:
                ax.set_xscale(xscale)
            if shaded:
                ax.axvspan(*shaded, **SHADE)
            ax.grid(alpha=0.2)
        if len(series) > 1:
            # The series legend rides in the ratio panel, whose lower left is
            # empty in every suite; the mean panel's is not.
            bottom.legend(loc="lower left")
        top.text(
            0.03,
            0.97,
            "\n".join(annotations),
            transform=top.transAxes,
            va="top",
            fontsize=ANNOTATION_FONTSIZE,
            bbox=ANNOTATION_BOX,
        )
    axes[0, 0].set_xticks([1, 5, 10, 15, 20])

    mean_limits = _symmetric_limits(np.concatenate(mean_content), 0.0)
    ratio_limits = _symmetric_limits(np.concatenate(ratio_content), 1.0)
    for column in range(3):
        axes[0, column].set_ylim(mean_limits)
        axes[1, column].set_ylim(ratio_limits)
    axes[0, 0].set_ylabel(r"$\langle B\rangle / \sigma^\mathrm{ana}$")
    axes[1, 0].set_ylabel(r"$\sigma_\mathrm{emp} / \sigma^\mathrm{ana}$")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stats",
        required=True,
        help="The campaign's per-realisation statistics npz.",
    )
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args()
    build_figure(read_statistics(args.stats), Path(args.out_file))


if __name__ == "__main__":
    main()
