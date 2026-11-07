"""COSEBI B_n / sigma_n across catalog versions, at the fiducial cut and full range.

Each version's modes come from the modes bundle and its sigma_n from the
covariance bundle, so the figure computes nothing. Statistical evidence (PTEs)
is in cosebis_pte_matrix.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from plotting_utils import (
    ERRORBAR_DEFAULTS,
    FIG_WIDTH_FULL,
    MARKER_STYLES,
    PAPER_MPLSTYLE,
    draw_normalized_boxes_linear_scale,
    find_fiducial_index,
    X_OFFSETS,
    get_version_alpha,
    version_label,
)
from print_width import savefig_at_width

# The paper stacks Fig. 8's three version-comparison panels at 0.88\textwidth.
STACKED_PRINT_WIDTH = 0.88 * FIG_WIDTH_FULL
from scale_cuts import add_scale_cut_convention_argument, resolve_fiducial_scale_cuts

from catalog_versions import add_catalog_arguments, catalogs_from
from decision_values import FIDUCIAL_BLIND, MODE_COUNT, REPORTING_BINNING

plt.style.use(PAPER_MPLSTYLE)


def _create_stacked_bmode_figure(
    fiducial_datasets,
    full_datasets,
    nmodes,
    scale_cuts,
    fiducial_idx,
    y_limits=None,
):
    """Create a vertically stacked B-mode COSEBIS comparison figure.

    Plots B_n / sigma_n (dimensionless, in units of standard deviation).
    Top panel: full range, bottom panel: fiducial scale cut.
    For each mode, a box spans the range of all versions' error bars,
    with a line marking the fiducial version's value.

    """
    x_offsets = np.array(X_OFFSETS)

    fig, axes = plt.subplots(
        2, 1, figsize=(FIG_WIDTH_FULL, FIG_WIDTH_FULL * 0.385), sharex=True
    )

    modes = np.arange(1, nmodes + 1)
    scale_labels = {"fiducial": "Fiducial", "full": "Full"}

    panels = [
        ("full", full_datasets, axes[0], scale_cuts["full"]),
        ("fiducial", fiducial_datasets, axes[1], scale_cuts["fiducial"]),
    ]

    legend_handles = []
    legend_labels = []

    for panel_idx, (scale_key, datasets, ax, scale_cut) in enumerate(panels):
        # Draw version spread boxes (before data points)
        draw_normalized_boxes_linear_scale(
            ax,
            modes,
            datasets,
            y_norm_key="Bn_normalized",
            fiducial_idx=fiducial_idx,
            x_offsets=x_offsets,
        )

        for i, data in enumerate(datasets):
            offset = x_offsets[i] if i < len(x_offsets) else 0
            marker = data.get(
                "marker", MARKER_STYLES[i] if i < len(MARKER_STYLES) else "o"
            )
            fillstyle = data.get("fillstyle", "full")
            mfc = data["color"] if fillstyle == "full" else "none"
            line = ax.errorbar(
                modes + offset,
                data["Bn_normalized"],
                yerr=np.ones_like(data["Bn_normalized"]),
                fmt=marker,
                color=data["color"],
                alpha=data["alpha"],
                markerfacecolor=mfc,
                markeredgecolor=data["color"],
                **ERRORBAR_DEFAULTS,
                zorder=2,
            )
            if panel_idx == 0:
                legend_handles.append(line)
                legend_labels.append(data["label"])

        ax.axhline(0, color="black", linewidth=0.8, alpha=0.6)
        ax.axvspan(0.5, 6.5, color="0.95", alpha=0.5, zorder=0)
        ax.set_ylabel(r"$B_n / \sigma_n$")
        ax.set_xlim(0.5, nmodes + 0.5)
        if y_limits:
            ax.set_ylim(y_limits)
        ax.tick_params(axis="both", width=0.5, length=3)

        label = scale_labels[scale_key]
        ax.text(
            0.98,
            0.95,
            rf"{label} $\theta = {scale_cut[0]:.0f}$--${scale_cut[1]:.0f}'$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    axes[1].set_xlabel("COSEBIS mode $n$")
    axes[1].set_xticks(np.arange(1, nmodes + 1))

    axes[0].legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        ncol=4,
        frameon=True,
        framealpha=0.9,
        handletextpad=0.3,
        columnspacing=1.0,
    )

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.08)
    return fig


def _modes_npz(results_dir, ver):
    return f"{results_dir}/{ver}_cosebis_modes.npz"


def _covariance_npz(cov_dir, ver, blind):
    return f"{cov_dir}/{ver}_{blind}_cosebis_covariance.npz"


def main(
    catalogs,
    blind,
    nmodes,
    reporting,
    results_dir,
    cov_dir,
    out_file,
    scale_cut_convention,
    fiducial_version=None,
    fiducial_modes_path=None,
    fiducial_covariance_path=None,
):
    # Leak-corrected, non-ecut versions
    versions = catalogs.comparison
    version_labels = catalogs.labels

    # Fiducial version whose inputs may be overridden with explicit paths
    fiducial_version = fiducial_version or catalogs.fiducial

    # Per-version modes and covariance, with explicit-path overrides for the
    # fiducial version.
    modes_paths = [
        fiducial_modes_path
        if v == fiducial_version and fiducial_modes_path
        else _modes_npz(results_dir, v)
        for v in versions
    ]
    covariance_paths = [
        fiducial_covariance_path
        if v == fiducial_version and fiducial_covariance_path
        else _covariance_npz(cov_dir, v, blind)
        for v in versions
    ]

    # Which version gets the fiducial reference line in boxes
    fiducial_for_comparison = catalogs.fiducial

    fiducial_scale_cut = resolve_fiducial_scale_cuts(reporting, scale_cut_convention)
    full_scale_cut = (reporting.min_sep, reporting.max_sep)

    scale_cuts = {
        "fiducial": fiducial_scale_cut,
        "full": full_scale_cut,
    }

    # Color/marker assignment: pair by parent version if ecut versions present
    has_ecut = any("_ecut" in v for v in versions)
    if has_ecut:
        parents = []
        for v in versions:
            parent = v.split("_ecut")[0].replace("_leak_corr", "")
            if parent not in parents:
                parents.append(parent)
        parent_colors = sns.husl_palette(len(parents), l=0.4)
        parent_color_map = dict(zip(parents, parent_colors))
    else:
        colors = sns.husl_palette(len(versions), l=0.4)

    all_datasets = {}

    for scale_key in scale_cuts:
        datasets = []

        for i, (version, modes_path, covariance_path) in enumerate(
            zip(versions, modes_paths, covariance_paths)
        ):
            if has_ecut:
                parent = version.split("_ecut")[0].replace("_leak_corr", "")
                color = parent_color_map[parent]
                is_cut = "_ecut" in version
                marker = "o" if not is_cut else "D"
                fillstyle = "full" if not is_cut else "none"
            else:
                color = colors[i]
                marker = MARKER_STYLES[i] if i < len(MARKER_STYLES) else "o"
                fillstyle = "full"

            Bn = np.load(modes_path, allow_pickle=True)[f"{scale_key}_Bn"]
            sigma_B = np.load(covariance_path, allow_pickle=True)[f"{scale_key}_sigma_B"]

            datasets.append(
                {
                    "version": version,
                    "label": version_label(version, version_labels),
                    "color": color,
                    "marker": marker,
                    "fillstyle": fillstyle,
                    "alpha": get_version_alpha(
                        version, fiducial_for_comparison
                    ),
                    "Bn_normalized": Bn / sigma_B,
                }
            )

        all_datasets[scale_key] = datasets

    # Compute y-limits
    all_y = []
    for datasets in all_datasets.values():
        for d in datasets:
            all_y.extend(d["Bn_normalized"] - 1)
            all_y.extend(d["Bn_normalized"] + 1)
    y_range = max(all_y) - min(all_y)
    y_limits = (min(all_y) - 0.1 * y_range, max(all_y) + 0.1 * y_range)

    # Find fiducial version index for box highlighting (use fiducial datasets)
    fiducial_idx = find_fiducial_index(
        all_datasets["fiducial"], fiducial_for_comparison
    )

    # Create output
    fig = _create_stacked_bmode_figure(
        all_datasets["fiducial"],
        all_datasets["full"],
        nmodes,
        scale_cuts,
        fiducial_idx,
        y_limits,
    )

    savefig_at_width(fig, out_file, STACKED_PRINT_WIDTH, dpi=300)
    print(f"Saved {out_file}")
    plt.close(fig)


def _from_cli(argv=None):

    ap = argparse.ArgumentParser(
        description="COSEBI B-mode version-comparison figure (4 leak-corrected catalogs)."
    )
    add_catalog_arguments(ap)
    ap.add_argument("--blind", type=FIDUCIAL_BLIND.__getitem__, required=True)
    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)
    ap.add_argument("--mode-count", type=MODE_COUNT.__getitem__, required=True)
    ap.add_argument(
        "--results-dir",
        required=True,
        help="Directory holding per-version <version>_cosebis_modes.npz files",
    )
    ap.add_argument(
        "--cov-dir",
        required=True,
        help="Directory holding per-version <version>_<blind>_cosebis_covariance.npz files",
    )
    ap.add_argument(
        "--out-file",
        required=True,
        help="Output file (lc {output}): the single stacked PNG",
    )
    ap.add_argument(
        "--fiducial-version",
        default=None,
        help="Version whose inputs may be overridden by --fiducial-modes-path "
        "and --fiducial-covariance-path (default: the fiducial version)",
    )
    ap.add_argument(
        "--fiducial-modes-path",
        default=None,
        help="The fiducial version's cosebis_modes_data npz",
    )
    ap.add_argument(
        "--fiducial-covariance-path",
        default=None,
        help="The fiducial version's cosebis_covariance npz",
    )
    add_scale_cut_convention_argument(ap)
    a = ap.parse_args(argv)
    main(
        catalogs_from(a),
        a.blind,
        a.mode_count,
        a.reporting_binning,
        a.results_dir,
        a.cov_dir,
        a.out_file,
        a.scale_cut_convention,
        fiducial_version=a.fiducial_version,
        fiducial_modes_path=a.fiducial_modes_path,
        fiducial_covariance_path=a.fiducial_covariance_path,
    )


if __name__ == "__main__":
    _from_cli()
