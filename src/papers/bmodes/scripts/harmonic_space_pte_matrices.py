"""Harmonic-space PTE matrix figures for B-modes paper.

Produces:
- Results: single-panel fiducial Cl^BB PTE matrix
- Appendix: N-panel composite for all versions from the catalog roster

Uses fiducial blind covariance only (blind independence validated elsewhere).
"""

import argparse
import json
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numerics import compute_chi2_pte, ell_bin_index, resolve_pte_threshold
from plotting_utils import (
    FIG_WIDTH_FULL,
    FIG_WIDTH_SINGLE,
    PAPER_MPLSTYLE,
    format_pte_colorbar,
    make_pte_colormap,
    make_pte_norm,
)
from print_width import savefig_at_width
from catalog_versions import VERSIONS, add_catalog_arguments, catalogs_from
from decision_values import FIDUCIAL_BLIND, FIDUCIAL_ELL_CUT

plt.style.use(PAPER_MPLSTYLE)


def _cl(results_dir, ver, blind="A", nbins=32):
    return f"{results_dir}/cl_{ver}_blind={blind}_powspace_nbins={nbins}.fits"


def _cl_cov(results_dir, ver, blind="A", nbins=32):
    return (
        f"{results_dir}/cl_cov_{ver}_blind={blind}_powspace_nbins={nbins}.fits"
    )


def compute_pte_matrix(
    cl_path, cl_cov_path, fiducial_ell_min=None, fiducial_ell_max=None
):
    """Compute PTE matrix for all ell-bin cut combinations.

    Parameters
    ----------
    cl_path : str
        Path to Cl FITS file.
    cl_cov_path : str
        Path to Cl covariance FITS file.
    fiducial_ell_min : float, optional
        Minimum ell for fiducial cut.
    fiducial_ell_max : float, optional
        Maximum ell for fiducial cut.

    Returns
    -------
    pte_matrix : ndarray
        n_ell x n_ell PTE matrix.
    ell : ndarray
        Ell bin centers.
    stats : dict
        Summary statistics.
    """
    # Load Cl data
    hdu = fits.open(cl_path)
    data = hdu["PSEUDO_CELL"].data
    hdu.close()

    ell = data["ELL"]
    cl_bb = data["BB"]
    n_ell = len(ell)

    # Load covariance
    hdu_cov = fits.open(cl_cov_path)
    cov_bb = hdu_cov["COVAR_BB_BB"].data
    hdu_cov.close()

    # Compute PTE matrix
    pte_matrix = np.full((n_ell, n_ell), np.nan)

    for i_min in range(n_ell):
        for i_max in range(i_min + 1, n_ell):
            idx_slice = slice(i_min, i_max + 1)
            bb_slice = cl_bb[idx_slice]
            cov_slice = cov_bb[idx_slice, idx_slice]

            try:
                chi2, pte, dof = compute_chi2_pte(bb_slice, cov_slice)
                pte_matrix[i_min, i_max] = pte
            except np.linalg.LinAlgError:
                pass

    # Compute statistics
    valid_ptes = pte_matrix[~np.isnan(pte_matrix)]
    pte_full_range = pte_matrix[0, n_ell - 1]

    stats_out = {
        "pte_at_full_range": float(pte_full_range),
        "n_evaluated": int(len(valid_ptes)),
        "ell_range": [float(ell.min()), float(ell.max())],
        "n_ell_bins": int(n_ell),
        "mean": float(np.nanmean(valid_ptes)),
        "median": float(np.nanmedian(valid_ptes)),
        "std": float(np.nanstd(valid_ptes)),
    }

    # Compute PTE at fiducial ell cuts if provided
    # Use bin edges: include only bins fully within [ell_min, ell_max]
    if fiducial_ell_min is not None and fiducial_ell_max is not None:
        i_fid_min, i_fid_max = ell_bin_index(ell, fiducial_ell_min, fiducial_ell_max)
        pte_fiducial = pte_matrix[i_fid_min, i_fid_max]
        stats_out["pte_at_fiducial"] = float(pte_fiducial)
        stats_out["fiducial_ell_range"] = [fiducial_ell_min, fiducial_ell_max]

    return pte_matrix, ell, stats_out


def plot_cl_pte_panel(
    ax,
    pte_matrix,
    ell,
    title,
    show_colorbar=False,
    show_xlabel=True,
    show_ylabel=True,
    fid_i_min=None,
    fid_i_max=None,
    pte_threshold=0.05,
):
    """Plot a single Cl PTE heatmap panel.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    pte_matrix : ndarray
        PTE matrix.
    ell : ndarray
        Ell bin centers.
    title : str
        Panel title.
    show_colorbar : bool
        Whether to add individual colorbar.
    show_xlabel, show_ylabel : bool
        Whether to show axis labels.

    Returns
    -------
    im : AxesImage
        The image object.
    """
    n_ell = len(ell)

    # Discrete colormap: solid blue below the selected threshold, solid red
    # above its complement, gradient between.
    pte_cmap = make_pte_colormap(pte_threshold)
    pte_norm = make_pte_norm(pte_threshold)

    # pte_matrix[i_min, i_max] -> transpose so rows=i_max, cols=i_min
    # With origin="lower", row 0 is at bottom (small i_max), row n-1 at top (large i_max)
    # This matches tick labels: y-axis has small ell at bottom, large ell at top
    pte_plot = pte_matrix.T

    im = ax.imshow(
        pte_plot,
        origin="lower",
        aspect="equal",
        cmap=pte_cmap,
        norm=pte_norm,
        extent=[0, n_ell, 0, n_ell],
    )

    # Mark fiducial ell cut (black square)
    if fid_i_min is not None and fid_i_max is not None:
        ax.add_patch(
            Rectangle(
                (fid_i_min, fid_i_max),
                1,
                1,
                fill=False,
                edgecolor="black",
                linewidth=1.5,
            )
        )

    # Tick positioning: x-axis at left edge of bins, y-axis at top edge
    tick_step = max(1, n_ell // 5)  # Fewer ticks for cleaner labels
    tick_indices = np.arange(0, n_ell, tick_step)

    # x-axis (lower cut): ticks at left edge of bins
    ax.set_xticks(tick_indices)
    # y-axis (upper cut): ticks at top edge of bins
    ax.set_yticks(tick_indices + 1)

    if show_xlabel:
        ax.set_xticklabels(
            [f"{ell[i]:.0f}" for i in tick_indices],
            fontsize=9,
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
    else:
        ax.set_xticklabels([])

    if show_ylabel:
        y_tick_labels = [f"{ell[min(i + 1, n_ell - 1)]:.0f}" for i in tick_indices]
        ax.set_yticklabels(y_tick_labels, fontsize=9)
    else:
        ax.set_yticklabels([])

    if show_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax, spacing="proportional")
        format_pte_colorbar(cbar, pte_threshold)
        cbar.set_label("PTE", fontsize=10)
        cbar.ax.tick_params(labelsize=9)

    return im


def create_single_panel(
    pte_matrix, ell, fid_i_min=None, fid_i_max=None, pte_threshold=0.05
):
    """Create single-panel figure for fiducial version."""
    fig, ax = plt.subplots(1, 1, figsize=(FIG_WIDTH_SINGLE, FIG_WIDTH_SINGLE))

    plot_cl_pte_panel(
        ax,
        pte_matrix,
        ell,
        "",
        show_colorbar=True,
        fid_i_min=fid_i_min,
        fid_i_max=fid_i_max,
        pte_threshold=pte_threshold,
    )

    ax.set_xlabel(r"$\ell_{\mathrm{min}}$")
    ax.set_ylabel(r"$\ell_{\mathrm{max}}$")

    plt.tight_layout()
    return fig


def create_npanel_composite(
    matrices, ells, panel_labels, fid_indices=None, pte_threshold=0.05
):
    """Create N-panel composite for appendix versions.

    Parameters
    ----------
    matrices : list of ndarray
        PTE matrices for each version.
    ells : list of ndarray
        Ell arrays for each version.
    panel_labels : list of str
        Labels for each panel.

    Returns
    -------
    fig : Figure
        The composite figure.
    """
    n_panels = len(matrices)
    gs_left, gs_right = 0.08, 0.94
    gs_bottom, gs_top = 0.28, 0.92
    wspace_val = 0.08
    cbar_ratio = 0.04
    fig_width = 7.0

    # Compute fig_height so panels are exactly square — no vertical whitespace
    total_width_units = n_panels + cbar_ratio + n_panels * wspace_val
    panel_frac = (gs_right - gs_left) / total_width_units
    panel_in = panel_frac * fig_width
    fig_height = panel_in / (gs_top - gs_bottom)

    fig = plt.figure(figsize=(fig_width, fig_height))
    width_ratios = [1] * n_panels + [cbar_ratio]
    gs = fig.add_gridspec(
        1,
        n_panels + 1,
        width_ratios=width_ratios,
        wspace=wspace_val,
        left=gs_left,
        right=gs_right,
        bottom=gs_bottom,
        top=gs_top,
    )

    axes = [fig.add_subplot(gs[i]) for i in range(n_panels)]
    cax_placeholder = fig.add_subplot(gs[n_panels])  # position reference

    for i, (ax, matrix, ell, label) in enumerate(
        zip(axes, matrices, ells, panel_labels)
    ):
        fi = fid_indices[i] if fid_indices else (None, None)
        im = plot_cl_pte_panel(
            ax,
            matrix,
            ell,
            "",
            show_ylabel=(i == 0),
            fid_i_min=fi[0],
            fid_i_max=fi[1],
            pte_threshold=pte_threshold,
        )
        # Version label inside panel (bottom-right)
        ax.text(
            0.95,
            0.05,
            label,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
        )

    # Draw first to get actual panel positions after aspect="equal" constraint
    fig.canvas.draw()
    ax0_pos = axes[0].get_position()
    cax_pos = cax_placeholder.get_position()
    fig.delaxes(cax_placeholder)

    # Colorbar at exact rendered panel height
    cax = fig.add_axes([cax_pos.x0, ax0_pos.y0, cax_pos.width, ax0_pos.height])
    cbar = fig.colorbar(im, cax=cax, spacing="proportional")
    format_pte_colorbar(cbar, pte_threshold)
    cbar.set_label("PTE", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Common axis labels
    fig.text(
        0.5 * (gs_left + gs_right),
        0.4 * gs_bottom,
        r"$\ell_{\mathrm{min}}$",
        ha="center",
        va="center",
    )
    fig.text(
        0.25 * gs_left,
        0.5 * (gs_bottom + gs_top),
        r"$\ell_{\mathrm{max}}$",
        va="center",
        ha="center",
        rotation="vertical",
    )

    return fig


def _matrix_to_json(matrix):
    """NaN-safe nested-list encoding of a PTE matrix (NaN -> null)."""
    return [[None if np.isnan(v) else float(v) for v in row] for row in matrix]


def _matrix_from_json(rows):
    """Inverse of _matrix_to_json (null -> NaN)."""
    return np.array(
        [[np.nan if v is None else float(v) for v in row] for row in rows], dtype=float
    )


def compute_all(
    catalogs,
    fiducial_blind,
    results_dir,
    fiducial_version_override=None,
    fiducial_cl_path=None,
    fiducial_cl_cov_path=None,
):
    """Run the per-version PTE sweep; return (matrices, ells, stats) dicts."""
    versions = [v for v in VERSIONS if "_ecut" not in v]

    fiducial_ell_min, fiducial_ell_max = FIDUCIAL_ELL_CUT["adopted"]

    # Per-version Cl / covariance files (canonical COSMO_VAL tree, blind A).
    # When --fiducial-version matches and an explicit override path is set, read
    # that path instead of the reconstructed pattern (override files lack the
    # blind=/powspace_nbins= tokens and are distinguished by directory).
    version_to_cl = {
        ver: (
            fiducial_cl_path
            if ver == fiducial_version_override and fiducial_cl_path is not None
            else _cl(results_dir, ver, blind=fiducial_blind)
        )
        for ver in versions
    }
    version_to_cov = {
        ver: (
            fiducial_cl_cov_path
            if ver == fiducial_version_override
            and fiducial_cl_cov_path is not None
            else _cl_cov(results_dir, ver, blind=fiducial_blind)
        )
        for ver in versions
    }

    # Compute PTE matrices for all versions using fiducial blind
    all_stats = {}
    all_matrices = {}
    all_ells = {}

    for version in versions:
        print(f"\n--- Processing {version} (blind={fiducial_blind}) ---")

        if version not in version_to_cl:
            print(f"  WARNING: No Cl file found for {version}, skipping")
            continue

        if version not in version_to_cov:
            print(f"  WARNING: No covariance file found for {version}, skipping")
            continue

        pte_matrix, ell, stats = compute_pte_matrix(
            version_to_cl[version],
            version_to_cov[version],
            fiducial_ell_min=fiducial_ell_min,
            fiducial_ell_max=fiducial_ell_max,
        )

        all_stats[version] = stats
        all_matrices[version] = pte_matrix
        all_ells[version] = ell

        print(f"  PTE at fiducial: {stats.get('pte_at_fiducial', 'N/A'):.4f}")
        print(f"  PTE at full range: {stats.get('pte_at_full_range', 'N/A'):.4f}")

    return all_matrices, all_ells, all_stats


def emit_matrices(catalogs, fiducial_blind, all_matrices, all_ells, out_file):
    """Write the per-version Cl^BB PTE grid NPZ to exactly `out_file`."""
    fiducial_version = catalogs.fiducial
    fiducial_ell_min, fiducial_ell_max = FIDUCIAL_ELL_CUT["adopted"]

    npz_payload = {
        "versions": np.array(list(all_matrices.keys())),
        "fiducial_version": fiducial_version,
        "fiducial_blind": fiducial_blind,
        "fiducial_ell_min": fiducial_ell_min,
        "fiducial_ell_max": fiducial_ell_max,
    }
    for version in all_matrices:
        npz_payload[f"{version}__pte_matrix"] = all_matrices[version]
        npz_payload[f"{version}__ell"] = all_ells[version]
    with open(out_file, "wb") as f:
        np.savez(f, **npz_payload)
    print(f"Saved PTE matrices to {out_file}")


def emit_evidence(catalogs, fiducial_blind, all_matrices, all_ells, all_stats, out_file):
    """Write the evidence JSON, carrying the FULL per-version PTE grids.

    The grids (and their ell axes) are in the payload so that the fiducial and
    appendix figures can be rendered from this file alone, recomputing nothing.
    """
    fiducial_version = catalogs.fiducial
    fiducial_ell_min, fiducial_ell_max = FIDUCIAL_ELL_CUT["adopted"]

    evidence_data = {
        "spec_id": "harmonic_space_pte_matrices",
        "generated": datetime.now().isoformat(),
        "evidence": {
            "blind": fiducial_blind,
            "fiducial_version": fiducial_version,
            "fiducial_ell_min": fiducial_ell_min,
            "fiducial_ell_max": fiducial_ell_max,
            "versions": {},
        },
    }

    for version, stats in all_stats.items():
        role = "fiducial" if version == fiducial_version else "appendix"
        evidence_data["evidence"]["versions"][version] = {
            "role": role,
            **stats,
            "ell": [float(v) for v in all_ells[version]],
            "pte_matrix": _matrix_to_json(all_matrices[version]),
        }

    with open(out_file, "w") as f:
        json.dump(evidence_data, f, indent=2)
    print(f"Saved evidence to {out_file}")


def _figure_payload_from_evidence(evidence_path):
    """Load (matrices, ells, fiducial metadata) back out of an evidence JSON."""
    with open(evidence_path) as f:
        payload = json.load(f)
    ev = payload["evidence"]
    versions = ev["versions"]
    missing = [v for v, d in versions.items() if "pte_matrix" not in d]
    if missing:
        raise SystemExit(
            "evidence JSON carries no per-version PTE grids for: "
            + ", ".join(missing)
            + " — regenerate it with --emit evidence"
        )
    all_matrices = {v: _matrix_from_json(d["pte_matrix"]) for v, d in versions.items()}
    all_ells = {v: np.array(d["ell"], dtype=float) for v, d in versions.items()}
    return all_matrices, all_ells, ev


def render_figure(
    catalogs, all_matrices, all_ells, which, out_file, meta=None, pte_threshold=0.05
):
    """Render the fiducial single panel or the appendix composite to `out_file`."""
    meta = meta or {}
    fiducial_version = meta.get("fiducial_version", catalogs.fiducial)
    fiducial_ell_min = meta.get("fiducial_ell_min", FIDUCIAL_ELL_CUT["adopted"][0])
    fiducial_ell_max = meta.get("fiducial_ell_max", FIDUCIAL_ELL_CUT["adopted"][1])
    version_labels = catalogs.labels

    # Compute fiducial ell indices per version using bin edges
    version_fid_indices = {
        version: ell_bin_index(ell, fiducial_ell_min, fiducial_ell_max)
        for version, ell in all_ells.items()
    }

    if which == "fiducial":
        if fiducial_version not in all_matrices:
            raise SystemExit(
                f"fiducial version {fiducial_version} absent from the PTE matrices"
            )
        fi_fid = version_fid_indices.get(fiducial_version, (None, None))
        fig = create_single_panel(
            all_matrices[fiducial_version],
            all_ells[fiducial_version],
            fid_i_min=fi_fid[0],
            fid_i_max=fi_fid[1],
            pte_threshold=pte_threshold,
        )
        fig.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        # Appendix N-panel composite (paper versions only, no ecut variants)
        appendix_versions = [
            v for v in VERSIONS if v in all_matrices and v in version_labels
        ]
        if len(appendix_versions) < 2:
            raise SystemExit(
                "appendix composite needs at least two labelled versions, got "
                f"{len(appendix_versions)}"
            )
        fig = create_npanel_composite(
            [all_matrices[v] for v in appendix_versions],
            [all_ells[v] for v in appendix_versions],
            [version_labels.get(v, v) for v in appendix_versions],
            fid_indices=[
                version_fid_indices.get(v, (None, None)) for v in appendix_versions
            ],
            pte_threshold=pte_threshold,
        )
        # The paper prints the appendix composite at 0.95\textwidth.
        savefig_at_width(
            fig, out_file, 0.95 * FIG_WIDTH_FULL, dpi=300, facecolor="white"
        )
        plt.close(fig)

    print(f"Saved {which} figure to {out_file}")


def _from_cli(argv=None):

    ap = argparse.ArgumentParser(
        description=(
            "Harmonic-space Cl^BB PTE matrices over the (ell_min, ell_max) "
            "bandpower-cut grid. Compute once with --emit "
            "{evidence,matrices,figure} --out-file; render the fiducial panel "
            "or appendix composite with --from-evidence --figure --out-file."
        )
    )
    add_catalog_arguments(ap)
    ap.add_argument("--blind", type=FIDUCIAL_BLIND.__getitem__)
    ap.add_argument(
        "--results-dir",
        default=None,
        help="Flat dir of per-version cl_* / cl_cov_* FITS "
        "(e.g. an untarred cl_sweep bundle); required unless --from-evidence",
    )
    ap.add_argument(
        "--emit",
        choices=["evidence", "matrices", "figure"],
        default=None,
        help="Which single artifact the compute pass writes to --out-file",
    )
    ap.add_argument(
        "--from-evidence",
        default=None,
        help="Render a figure from an evidence JSON written by --emit evidence "
        "(recomputes nothing)",
    )
    ap.add_argument(
        "--figure",
        choices=["fiducial", "appendix"],
        default="fiducial",
        help="Which figure to render (--emit figure or --from-evidence)",
    )
    ap.add_argument(
        "--out-file", required=True, help="Path of the single output file (lc {output})"
    )
    ap.add_argument(
        "--pte-threshold",
        type=resolve_pte_threshold,
        default="nominal_0p05",
        help="ASTRA PTE-threshold option ID or numeric threshold",
    )
    ap.add_argument(
        "--fiducial-version",
        default=None,
        help="Version whose inputs are read from explicit lc override paths below",
    )
    ap.add_argument(
        "--fiducial-cl-path",
        default=None,
        help="Explicit path to fiducial Cl FITS produced by lc "
        "(overrides pattern reconstruction for --fiducial-version)",
    )
    ap.add_argument(
        "--fiducial-cl-cov-path",
        default=None,
        help="Explicit path to fiducial Cl covariance FITS produced by lc "
        "(overrides pattern reconstruction for --fiducial-version)",
    )
    a = ap.parse_args(argv)
    catalogs = catalogs_from(a)

    if a.from_evidence:
        all_matrices, all_ells, meta = _figure_payload_from_evidence(a.from_evidence)
        render_figure(
            catalogs,
            all_matrices,
            all_ells,
            a.figure,
            a.out_file,
            meta=meta,
            pte_threshold=a.pte_threshold,
        )
        return

    if a.emit is None:
        ap.error("one of --emit or --from-evidence is required")
    if a.results_dir is None:
        ap.error("--results-dir is required for --emit")

    all_matrices, all_ells, all_stats = compute_all(
        catalogs,
        a.blind,
        a.results_dir,
        fiducial_version_override=a.fiducial_version,
        fiducial_cl_path=a.fiducial_cl_path,
        fiducial_cl_cov_path=a.fiducial_cl_cov_path,
    )

    if a.emit == "matrices":
        emit_matrices(catalogs, a.blind, all_matrices, all_ells, a.out_file)
    elif a.emit == "evidence":
        emit_evidence(catalogs, a.blind, all_matrices, all_ells, all_stats, a.out_file)
    else:
        render_figure(
            catalogs,
            all_matrices,
            all_ells,
            a.figure,
            a.out_file,
            pte_threshold=a.pte_threshold,
        )


if __name__ == "__main__":
    _from_cli()
