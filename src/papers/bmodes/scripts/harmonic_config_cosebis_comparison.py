"""Harmonic vs configuration-space COSEBIS cross-validation.

Cross-validates COSEBIS E_n and B_n computed from two independent paths:
1. Harmonic: pseudo-C_ell (powspace) -> cosebis_from_Cell() with 100k-point FFT-log
2. Configuration: fine-binned xi_pm -> calculate_cosebis()

Produces 9 figures per the standard data vector pattern (1 paper + 4 corrected
+ 4 uncorrected), plus a version comparison figure.

Angular range parameterized by ``--angular-range`` (full or fiducial).
Modes > reliable_mode_max shown with gray band.
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import seaborn as sns
import treecorr
from astropy.io import fits
from cosmo_numba.B_modes.cosebis import COSEBIS
from numerics import compute_chi2_pte
from plotting_utils import (
    EXCLUDED_SHADE,
    FIG_WIDTH_FULL,
    FIG_WIDTH_SINGLE,
    PAPER_MPLSTYLE,
    iter_version_figures,
    version_label,
)
from scale_cuts import add_scale_cut_convention_argument, resolve_fiducial_scale_cuts

from sp_validation.b_modes import calculate_cosebis
from catalog_versions import VERSIONS, catalogs_from
from decision_values import (
    FIDUCIAL_BLIND,
    FIDUCIAL_CATALOG,
    HARMONIC_VALIDATION_BINNING,
    INTEGRATION_GRID,
    MODE_COUNT,
    PATCH_COUNT,
    REPORTING_BINNING,
)

plt.style.use(PAPER_MPLSTYLE)

# At 32 bins, modes > 6 are unreliable. At 96 bins, modes 1-6 are recovered
# to <0.3% on GLASS mocks. Higher modes are progressively less reliable due
# to W_n(ell) oscillation resolution. We use 6 to match the mode range used
# for cosmological inference throughout the paper.
_RELIABLE_MODE_MAX_BY_NBINS = {32: 6, 96: 6, 128: 6}


def _safe_sigma(covariance_diag):
    """Standard deviation from covariance diagonal, clamped to avoid division by zero."""
    sigma = np.sqrt(np.clip(covariance_diag, 0, None))
    return np.where(sigma > 0, sigma, 1.0)


def _build_cell_covariance(cov_path):
    """Build joint [EE, BB] covariance from FITS pseudo-Cl covariance file."""
    with fits.open(cov_path) as hdul:
        cov_EE_EE = hdul["COVAR_EE_EE"].data
        cov_EE_BB = hdul["COVAR_EE_BB"].data
        cov_BB_EE = hdul["COVAR_BB_EE"].data
        cov_BB_BB = hdul["COVAR_BB_BB"].data

    n_ell = cov_EE_EE.shape[0]
    cov = np.zeros((2 * n_ell, 2 * n_ell))
    cov[:n_ell, :n_ell] = cov_EE_EE
    cov[:n_ell, n_ell:] = cov_EE_BB
    cov[n_ell:, :n_ell] = cov_BB_EE
    cov[n_ell:, n_ell:] = cov_BB_BB
    return cov


def _build_transforms(cosebis_obj, ell):
    """Build linear transform matrices T_E, T_B via basis-vector approach.

    T_E[n, i] = E_n response to unit C_ell in bin i (B=0).
    T_B[n, i] = B_n response to unit C_ell in bin i (E=0).
    """
    n_ell = len(ell)
    nmodes = cosebis_obj.Nmax
    zeros = np.zeros(n_ell)

    T_E = np.zeros((nmodes, n_ell))
    T_B = np.zeros((nmodes, n_ell))

    for idx in range(n_ell):
        basis = np.zeros(n_ell)
        basis[idx] = 1.0

        ce_basis, _ = cosebis_obj.cosebis_from_Cell(
            ell=ell,
            Cell_E=basis,
            Cell_B=zeros,
            cache=True,
        )
        _, cb_basis = cosebis_obj.cosebis_from_Cell(
            ell=ell,
            Cell_E=zeros,
            Cell_B=basis,
            cache=True,
        )

        T_E[:, idx] = ce_basis
        T_B[:, idx] = cb_basis

    return T_E, T_B


def _compute_harmonic_cosebis(
    pseudo_cl_path, pseudo_cov_path, nmodes, theta_min, theta_max
):
    """Compute COSEBIS from pseudo-C_ell and propagate covariance."""
    with fits.open(pseudo_cl_path) as hdul:
        data = hdul["PSEUDO_CELL"].data
        ell = np.asarray(data["ELL"], dtype=float)
        cl_ee = np.asarray(data["EE"], dtype=float)
        cl_bb = np.asarray(data["BB"], dtype=float)

    cosebis_obj = COSEBIS(theta_min, theta_max, nmodes)

    ce_harm, cb_harm = cosebis_obj.cosebis_from_Cell(
        ell=ell,
        Cell_E=cl_ee,
        Cell_B=cl_bb,
        cache=True,
    )

    T_E, T_B = _build_transforms(cosebis_obj, ell)

    cov_cell = _build_cell_covariance(pseudo_cov_path)
    n_ell = len(ell)
    transform = np.zeros((2 * nmodes, 2 * n_ell))
    transform[:nmodes, :n_ell] = T_E
    transform[nmodes:, n_ell:] = T_B
    cov_harmonic = transform @ cov_cell @ transform.T

    return ce_harm, cb_harm, cov_harmonic


def _compute_config_cosebis(xi_path, cov_path, nmodes, scale_cut, integration):
    """Compute COSEBIS from fine-binned xi_pm for a single scale cut.

    Returns (En, Bn, cov).
    """
    min_sep_int, max_sep_int, nbins_int = integration

    gg = treecorr.GGCorrelation(
        min_sep=min_sep_int,
        max_sep=max_sep_int,
        nbins=nbins_int,
        sep_units="arcmin",
    )
    gg.read(xi_path)

    results = calculate_cosebis(
        gg, nmodes=nmodes, scale_cuts=[scale_cut], cov_path=cov_path
    )
    r = results[scale_cut]
    return r["En"], r["Bn"], r["cov"]


def _shade_reliable(ax, nmodes, reliable_mode_max=6):
    """Shade the modes beyond reliable_mode_max, as the paper's other COSEBI figures do."""
    ax.axvspan(reliable_mode_max + 0.5, nmodes + 0.5, **EXCLUDED_SHADE)


def _make_data_vector_figure(
    harm_data, config_data, nmodes, scale_cut, title=None, reliable_mode_max=6
):
    """Data vector figure: E and B from both methods (single scale cut).

    Two rows (E-modes, B-modes).
    E-modes in raw units. B-modes in B_n/sigma units.
    Gray band over the modes beyond reliable_mode_max.
    """
    modes = np.arange(1, nmodes + 1)
    c_cfg, c_harm = sns.color_palette("Dark2", 2)

    ce_config, cb_config, cov_config = config_data
    sigma_config_E = np.sqrt(np.clip(np.diag(cov_config[:nmodes, :nmodes]), 0, None))
    sigma_config_B = _safe_sigma(np.diag(cov_config[nmodes:, nmodes:]))

    ce_harm, cb_harm = harm_data["En"], harm_data["Bn"]
    sigma_harm_E = harm_data["sigma_E"]
    sigma_harm_B = np.where(harm_data["sigma_B"] > 0, harm_data["sigma_B"], 1.0)

    fig, (ax_e, ax_b) = plt.subplots(
        2, 1, figsize=(FIG_WIDTH_SINGLE, FIG_WIDTH_SINGLE * 0.85), sharex=True
    )

    # --- E-modes (factored by 10^10 for readability) ---
    e_scale = 1e10
    ax_e.errorbar(
        modes - 0.1,
        ce_config * e_scale,
        yerr=sigma_config_E * e_scale,
        fmt="o",
        color=c_cfg,
        mfc=c_cfg,
        ms=4,
        alpha=0.8,
        capsize=2,
        capthick=0.8,
        elinewidth=0.8,
        label="configuration space",
    )
    ax_e.errorbar(
        modes + 0.1,
        ce_harm * e_scale,
        yerr=sigma_harm_E * e_scale,
        fmt="s",
        color=c_harm,
        mfc="white",
        mew=0.8,
        ms=4,
        alpha=0.8,
        capsize=2,
        capthick=0.8,
        elinewidth=0.8,
        label="harmonic space",
    )

    ax_e.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    _shade_reliable(ax_e, nmodes, reliable_mode_max)
    ax_e.set_ylabel(r"$E_n \times 10^{10}$")
    ax_e.set_xlim(0.5, nmodes + 0.5)
    ax_e.tick_params(axis="both", width=0.5, length=3)
    if title:
        ax_e.set_title(title)
    ax_e.legend(loc="upper center", framealpha=0.9)

    # --- B-modes (B_n / sigma units) ---
    ax_b.errorbar(
        modes - 0.1,
        cb_config / sigma_config_B,
        yerr=np.ones(nmodes),
        fmt="o",
        color=c_cfg,
        mfc=c_cfg,
        ms=4,
        alpha=0.8,
        capsize=2,
        capthick=0.8,
        elinewidth=0.8,
    )
    ax_b.errorbar(
        modes + 0.1,
        cb_harm / sigma_harm_B,
        yerr=np.ones(nmodes),
        fmt="s",
        color=c_harm,
        mfc="white",
        mew=0.8,
        ms=4,
        alpha=0.8,
        capsize=2,
        capthick=0.8,
        elinewidth=0.8,
    )

    ax_b.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    _shade_reliable(ax_b, nmodes, reliable_mode_max)
    ax_b.set_ylabel(r"$B_n / \sigma_n$")
    ax_b.set_xlabel("COSEBI mode $n$")
    ax_b.set_xticks(np.arange(1, nmodes + 1))
    ax_b.set_xlim(0.5, nmodes + 0.5)
    ax_b.tick_params(axis="both", width=0.5, length=3)

    ax_b.text(
        0.98,
        0.95,
        rf"$\theta = {scale_cut[0]:.0f}$--${scale_cut[1]:.0f}'$",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.15)
    return fig


def _make_combined_data_vector_figure(
    harm_data_full,
    config_data_full,
    harm_data_fid,
    config_data_fid,
    nmodes,
    scale_cut_full,
    scale_cut_fid,
    title=None,
    reliable_mode_max=6,
):
    """Combined 2x2 data vector figure: rows=E/B, columns=full/fiducial.

    Column titles show the theta range. Gray band over the modes beyond 6.
    Dark2 palette, filled circles for config-space, open squares for harmonic.
    """
    modes = np.arange(1, nmodes + 1)
    c_cfg = "#2c5f8a"  # steel blue (config-space)
    c_harm = "#c45a2c"  # burnt orange (harmonic-space)

    fig, axes = plt.subplots(
        2, 2, figsize=(FIG_WIDTH_FULL, FIG_WIDTH_FULL * 0.4), sharex=True
    )

    datasets = [
        (config_data_full, harm_data_full, scale_cut_full),
        (config_data_fid, harm_data_fid, scale_cut_fid),
    ]

    e_scale = 1e10

    for col, (cfg_data, h_data, sc) in enumerate(datasets):
        ce_cfg, cb_cfg, cov_cfg = cfg_data
        sigma_cfg_E = np.sqrt(np.clip(np.diag(cov_cfg[:nmodes, :nmodes]), 0, None))
        sigma_cfg_B = _safe_sigma(np.diag(cov_cfg[nmodes:, nmodes:]))
        ce_h, cb_h = h_data["En"], h_data["Bn"]
        sigma_h_E = h_data["sigma_E"]
        sigma_h_B = np.where(h_data["sigma_B"] > 0, h_data["sigma_B"], 1.0)

        # Column title
        col_titles = ["All scales", "Fiducial scale cuts"]
        axes[0, col].set_title(col_titles[col])

        # E-modes
        ax_e = axes[0, col]
        ax_e.errorbar(
            modes - 0.1,
            ce_cfg * e_scale,
            yerr=sigma_cfg_E * e_scale,
            fmt="o",
            color=c_cfg,
            mfc=c_cfg,
            ms=4,
            alpha=0.8,
            capsize=2,
            capthick=0.8,
            elinewidth=0.8,
            label="configuration space" if col == 0 else None,
        )
        ax_e.errorbar(
            modes + 0.1,
            ce_h * e_scale,
            yerr=sigma_h_E * e_scale,
            fmt="s",
            color=c_harm,
            mfc="white",
            mew=0.8,
            ms=4,
            alpha=0.8,
            capsize=2,
            capthick=0.8,
            elinewidth=0.8,
            label="harmonic space" if col == 0 else None,
        )
        ax_e.axhline(0.0, color="black", lw=0.8, alpha=0.6)
        _shade_reliable(ax_e, nmodes, reliable_mode_max)
        ax_e.set_xlim(0.5, nmodes + 0.5)
        ax_e.tick_params(axis="both", width=0.5, length=3)

        # B-modes
        ax_b = axes[1, col]
        ax_b.errorbar(
            modes - 0.1,
            cb_cfg / sigma_cfg_B,
            yerr=np.ones(nmodes),
            fmt="o",
            color=c_cfg,
            mfc=c_cfg,
            ms=4,
            alpha=0.8,
            capsize=2,
            capthick=0.8,
            elinewidth=0.8,
        )
        ax_b.errorbar(
            modes + 0.1,
            cb_h / sigma_h_B,
            yerr=np.ones(nmodes),
            fmt="s",
            color=c_harm,
            mfc="white",
            mew=0.8,
            ms=4,
            alpha=0.8,
            capsize=2,
            capthick=0.8,
            elinewidth=0.8,
        )
        ax_b.axhline(0.0, color="black", lw=0.8, alpha=0.6)
        _shade_reliable(ax_b, nmodes, reliable_mode_max)
        ax_b.set_xlabel("COSEBI mode $n$")
        ax_b.set_xticks(np.arange(1, nmodes + 1))
        ax_b.set_xticklabels(
            [str(i) for i in range(1, nmodes + 1)],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        ax_b.set_xlim(0.5, nmodes + 0.5)
        ax_b.tick_params(axis="both", width=0.5, length=3)

    # Fix B_n/sigma y-range (same for both columns, independent of E_n)
    for col in range(2):
        axes[1, col].set_ylim(-3.5, 5.5)
        axes[1, col].set_yticks([0, 3])

    # Row labels on left column only
    axes[0, 0].set_ylabel(r"$E_n \times 10^{10}$")
    axes[1, 0].set_ylabel(r"$B_n / \sigma_n$")

    # Legend in the top-left panel, centred over the shaded modes (where the
    # points sit near zero): x in data coordinates, y in axes coordinates.
    shaded_centre = 0.5 * ((reliable_mode_max + 0.5) + (nmodes + 0.5))
    axes[0, 0].legend(
        loc="upper center",
        bbox_to_anchor=(shaded_centre, 1.0),
        bbox_transform=mtransforms.blended_transform_factory(
            axes[0, 0].transData, axes[0, 0].transAxes
        ),
        framealpha=0.9,
    )

    if title:
        fig.suptitle(title)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.15, wspace=0.13)
    return fig


def _make_version_comparison_figure(
    all_version_data, nmodes, scale_cut, version_labels_map, reliable_mode_max=6
):
    """Version comparison: all versions, both methods, B-modes in B_n/sigma.

    Single panel, full angular range only.
    Config=filled marker, harmonic=open marker. Gray band over unreliable modes.
    """
    modes = np.arange(1, nmodes + 1)
    n_versions = len(all_version_data)
    colors = sns.color_palette("Dark2", n_versions)

    fig, ax = plt.subplots(figsize=(FIG_WIDTH_FULL, FIG_WIDTH_FULL * 0.35))
    x_offsets = np.linspace(-0.2, 0.2, n_versions)

    for i, (version, vdata) in enumerate(all_version_data.items()):
        color = colors[i]
        offset = x_offsets[i]
        label = version_label(version, version_labels_map)

        _, cb_config, cov_config = vdata["config"]
        sigma_config_B = _safe_sigma(np.diag(cov_config[nmodes:, nmodes:]))

        cb_harm = vdata["harm"]["Bn"]
        sigma_harm_B = np.where(
            vdata["harm"]["sigma_B"] > 0, vdata["harm"]["sigma_B"], 1.0
        )

        ax.errorbar(
            modes + offset - 0.03,
            cb_config / sigma_config_B,
            yerr=np.ones(nmodes),
            fmt="o",
            color=color,
            mfc=color,
            ms=3.5,
            alpha=0.85,
            capsize=1.5,
            capthick=0.6,
            elinewidth=0.6,
            label=f"{label} (config)",
        )
        ax.errorbar(
            modes + offset + 0.03,
            cb_harm / sigma_harm_B,
            yerr=np.ones(nmodes),
            fmt="s",
            color=color,
            mfc="white",
            ms=3.5,
            alpha=0.85,
            capsize=1.5,
            capthick=0.6,
            elinewidth=0.6,
            mew=0.8,
            label=f"{label} (harmonic)",
        )

    ax.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    _shade_reliable(ax, nmodes, reliable_mode_max)
    ax.set_ylabel(r"$B_n / \sigma_n$")
    ax.set_xlabel("COSEBI mode $n$")
    ax.set_xticks(np.arange(1, nmodes + 1))
    ax.set_xlim(0.5, nmodes + 0.5)
    ax.tick_params(axis="both", width=0.5, length=3)

    ax.text(
        0.98,
        0.95,
        rf"$\theta = {scale_cut[0]:.0f}$--${scale_cut[1]:.0f}'$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    ax.legend(
        loc="upper center",
        ncol=4,
        frameon=True,
        framealpha=0.9,
        handletextpad=0.3,
        columnspacing=0.8,
    )

    plt.tight_layout()
    return fig


def main(
    catalogs,
    nmodes,
    cosebis_nbins,
    integration,
    reporting,
    inputs,
    scale_cut,
    output_dir,
    scale_cut_convention,
    paper_figure_name=None,
    spec_path=None,
):
    """Harmonic-vs-config COSEBI cross-check for one angular range.

    ``inputs`` is a dict of per-version keys
    ``pseudo_cl_{ver}`` / ``pseudo_cl_cov_{ver}`` / ``xi_{ver}`` / ``cov_{ver}``
    (all absolute paths) plus ``specs``. All artifacts land under ``output_dir``;
    ``paper_figure_name`` (when set) is the combined 2×2 paper PDF filename.
    """
    version_labels_map = catalogs.labels

    fiducial_version = catalogs.fiducial
    scale_cut = tuple(scale_cut)

    reliable_mode_max = _reliable_mode_max(cosebis_nbins)
    print(f"Using {cosebis_nbins}-bin pseudo-Cl, reliable modes 1-{reliable_mode_max}")

    versions_leak_corr = [v for v in VERSIONS if "_leak_corr" in v]

    # Build input path lookups from the inputs dict
    pseudo_cl_paths = {
        k: v
        for k, v in inputs.items()
        if k.startswith("pseudo_cl_") and not k.startswith("pseudo_cl_cov_")
    }
    pseudo_cov_paths = {
        k: v for k, v in inputs.items() if k.startswith("pseudo_cl_cov_")
    }
    xi_paths = {k: v for k, v in inputs.items() if k.startswith("xi_")}
    cov_paths = {k: v for k, v in inputs.items() if k.startswith("cov_")}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Track generated artifacts
    output = {}

    # Compute B-mode PTEs for leak-corrected versions only
    harmonic_ptes = {}
    config_ptes = {}

    # Cache leak-corrected results for version comparison figure
    all_version_data = {}

    # Generate all 9 per-version data vector figures
    for fig_spec in iter_version_figures(version_labels_map, fiducial_version):
        if fig_spec["leak_corrected"]:
            ver = fig_spec["version_leak_corr"]
        else:
            ver = fig_spec["version_uncorr"]

        pseudo_cl_key = f"pseudo_cl_{ver}"
        pseudo_cov_key = f"pseudo_cl_cov_{ver}"
        xi_key = f"xi_{ver}"
        cov_key = f"cov_{ver}"

        print(
            f"Processing {ver} ({'corrected' if fig_spec['leak_corrected'] else 'uncorrected'})..."
        )

        # Config-space
        ce_cfg, cb_cfg, cov_cfg = _compute_config_cosebis(
            xi_paths[xi_key],
            cov_paths[cov_key],
            nmodes,
            scale_cut,
            integration,
        )

        # Harmonic-space
        ce_h, cb_h, cov_h = _compute_harmonic_cosebis(
            pseudo_cl_paths[pseudo_cl_key],
            pseudo_cov_paths[pseudo_cov_key],
            nmodes,
            scale_cut[0],
            scale_cut[1],
        )
        harm = {
            "En": ce_h,
            "Bn": cb_h,
            "sigma_E": np.sqrt(np.clip(np.diag(cov_h[:nmodes, :nmodes]), 0, None)),
            "sigma_B": np.sqrt(np.clip(np.diag(cov_h[nmodes:, nmodes:]), 0, None)),
        }

        # Create figure
        fig = _make_data_vector_figure(
            harm,
            (ce_cfg, cb_cfg, cov_cfg),
            nmodes,
            scale_cut,
            title=fig_spec["title"],
            reliable_mode_max=reliable_mode_max,
        )
        fig_path = output_dir / fig_spec["filename"]
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        print(f"  Saved {fig_path}")
        plt.close(fig)

        output[fig_spec["filename"].replace(".png", "")] = fig_spec["filename"]

        # For the paper figure, produce a combined 2x2 (both scale cuts)
        if fig_spec["is_paper_figure"] and paper_figure_name:
            # Compute the other scale cut for the combined figure
            scale_cut_full, scale_cut_fid = _scale_cuts(
                reporting, scale_cut_convention
            )
            # Determine which is "this" and which is "other"
            if scale_cut == scale_cut_full:
                other_cut = scale_cut_fid
            else:
                other_cut = scale_cut_full

            # Compute COSEBIS at the other scale cut for this version
            ce_cfg_other, cb_cfg_other, cov_cfg_other = _compute_config_cosebis(
                xi_paths[xi_key],
                cov_paths[cov_key],
                nmodes,
                other_cut,
                integration,
            )
            ce_h_other, cb_h_other, cov_h_other = _compute_harmonic_cosebis(
                pseudo_cl_paths[pseudo_cl_key],
                pseudo_cov_paths[pseudo_cov_key],
                nmodes,
                other_cut[0],
                other_cut[1],
            )
            harm_other = {
                "En": ce_h_other,
                "Bn": cb_h_other,
                "sigma_E": np.sqrt(
                    np.clip(np.diag(cov_h_other[:nmodes, :nmodes]), 0, None)
                ),
                "sigma_B": np.sqrt(
                    np.clip(np.diag(cov_h_other[nmodes:, nmodes:]), 0, None)
                ),
            }

            # Arrange: full on left, fiducial on right
            if scale_cut == scale_cut_full:
                harm_full, cfg_full = harm, (ce_cfg, cb_cfg, cov_cfg)
                harm_fid, cfg_fid = (
                    harm_other,
                    (ce_cfg_other, cb_cfg_other, cov_cfg_other),
                )
            else:
                harm_fid, cfg_fid = harm, (ce_cfg, cb_cfg, cov_cfg)
                harm_full, cfg_full = (
                    harm_other,
                    (ce_cfg_other, cb_cfg_other, cov_cfg_other),
                )

            fig_combined = _make_combined_data_vector_figure(
                harm_full,
                cfg_full,
                harm_fid,
                cfg_fid,
                nmodes,
                scale_cut_full,
                scale_cut_fid,
                title=fig_spec["title"] if fig_spec["title"] else None,
                reliable_mode_max=reliable_mode_max,
            )
            paper_path = output_dir / paper_figure_name
            paper_path.parent.mkdir(parents=True, exist_ok=True)
            fig_combined.savefig(paper_path, bbox_inches="tight")
            print(f"  Saved combined paper figure {paper_path}")
            plt.close(fig_combined)

            # cosebis_harmonic_modes data artifact (ASTRA cosebis.cosebis_harmonic_modes):
            # persist the fiducial-version harmonic + config-space COSEBI vectors and
            # covariances at both scale cuts.
            if scale_cut == scale_cut_full:
                cov_h_full, cov_h_fid = cov_h, cov_h_other
            else:
                cov_h_fid, cov_h_full = cov_h, cov_h_other
            np.savez(
                output_dir / f"cosebis_harmonic_modes_{ver}.npz",
                nmodes=nmodes,
                version=ver,
                cosebis_nbins=cosebis_nbins,
                reliable_mode_max=reliable_mode_max,
                full_scale_cut=np.array(scale_cut_full),
                fiducial_scale_cut=np.array(scale_cut_fid),
                full_harm_En=harm_full["En"],
                full_harm_Bn=harm_full["Bn"],
                full_harm_cov=cov_h_full,
                full_harm_sigma_E=harm_full["sigma_E"],
                full_harm_sigma_B=harm_full["sigma_B"],
                fiducial_harm_En=harm_fid["En"],
                fiducial_harm_Bn=harm_fid["Bn"],
                fiducial_harm_cov=cov_h_fid,
                fiducial_harm_sigma_E=harm_fid["sigma_E"],
                fiducial_harm_sigma_B=harm_fid["sigma_B"],
                full_config_En=cfg_full[0],
                full_config_Bn=cfg_full[1],
                full_config_cov=cfg_full[2],
                fiducial_config_En=cfg_fid[0],
                fiducial_config_Bn=cfg_fid[1],
                fiducial_config_cov=cfg_fid[2],
            )
            print(f"  Saved cosebis_harmonic_modes NPZ for {ver}")

        # Compute B-mode PTEs for leak-corrected versions
        if fig_spec["leak_corrected"] and ver not in harmonic_ptes:
            n_rel = reliable_mode_max
            # Harmonic-space PTE
            cb_h_rel = cb_h[:n_rel]
            cov_h_B_rel = cov_h[nmodes : nmodes + n_rel, nmodes : nmodes + n_rel]
            chi2, pte, dof = compute_chi2_pte(cb_h_rel, cov_h_B_rel)
            harmonic_ptes[ver] = {"chi2": chi2, "pte": pte, "dof": dof}
            print(
                f"  Harmonic B-mode PTE (modes 1-{n_rel}): {pte:.4f} (chi2={chi2:.2f}, dof={dof})"
            )
            # Config-space PTE
            cb_cfg_rel = cb_cfg[:n_rel]
            cov_cfg_B_rel = cov_cfg[nmodes : nmodes + n_rel, nmodes : nmodes + n_rel]
            chi2_c, pte_c, dof_c = compute_chi2_pte(cb_cfg_rel, cov_cfg_B_rel)
            config_ptes[ver] = {"chi2": chi2_c, "pte": pte_c, "dof": dof_c}
            print(
                f"  Config B-mode PTE (modes 1-{n_rel}): {pte_c:.4f} (chi2={chi2_c:.2f}, dof={dof_c})"
            )

        # Cache leak-corrected results for version comparison figure
        if fig_spec["leak_corrected"] and ver not in all_version_data:
            all_version_data[ver] = {
                "config": (ce_cfg, cb_cfg, cov_cfg),
                "harm": harm,
            }

    fig_vc = _make_version_comparison_figure(
        all_version_data,
        nmodes,
        scale_cut,
        version_labels_map,
        reliable_mode_max=reliable_mode_max,
    )
    fig_vc_path = output_dir / "figure_versions.png"
    fig_vc.savefig(fig_vc_path, dpi=300, bbox_inches="tight")
    print(f"Saved {fig_vc_path}")
    plt.close(fig_vc)
    output["figure_versions"] = "figure_versions.png"

    # --- Evidence ---
    evidence = {
        "spec_id": "harmonic_config_cosebis_comparison",
        "spec_path": spec_path
        or "papers/bmodes/config/harmonic_config_cosebis_comparison.md",
        "generated": datetime.now().isoformat(),
        "evidence": {
            "nmodes": nmodes,
            "reliable_modes": reliable_mode_max,
            "scale_cut": list(scale_cut),
            "powspace_nbins": cosebis_nbins,
            "versions_compared": versions_leak_corr,
            "harmonic_b_mode_ptes": harmonic_ptes,
            "config_b_mode_ptes": config_ptes,
            "note": (
                f"Harmonic-space COSEBIS computed via cosmo_numba cosebis_from_Cell() "
                f"from {cosebis_nbins}-bin powspace pseudo-C_ell. "
                f"With {cosebis_nbins} bins, modes 1-{reliable_mode_max} are reliably "
                f"recovered (validated on GLASS mocks: E_n/config within 2%). "
                f"Higher modes are limited by W_n(ell) numerical precision at "
                f"~10^-14 amplitudes, not binning. "
                f"Harmonic B-mode PTEs use modes 1-{reliable_mode_max} with "
                f"propagated covariance."
            ),
        },
        "output": output,
    }

    evidence_path = output_dir / "evidence.json"
    with open(evidence_path, "w") as f:
        json.dump(evidence, f, indent=2, default=str)
    print(f"Saved evidence to {evidence_path}")


def _scale_cuts(reporting, scale_cut_convention):
    """(full, fiducial) scale cuts, as used for the combined paper figure."""
    full = (reporting.min_sep, reporting.max_sep)
    fiducial = resolve_fiducial_scale_cuts(reporting, scale_cut_convention)
    return full, fiducial


def _reliable_mode_max(cosebis_nbins):
    return _RELIABLE_MODE_MAX_BY_NBINS.get(
        cosebis_nbins,
        8 if cosebis_nbins >= 64 else 6,
    )


def _compute_modes_npz(
    nmodes,
    cosebis_nbins,
    integration,
    reporting,
    version,
    pseudo_cl_path,
    pseudo_cov_path,
    xi_path,
    cov_path,
    out_file,
    scale_cut_convention,
):
    """Compute the harmonic- and config-space COSEBIs at both scale cuts for one
    version and write exactly one NPZ to ``out_file``."""
    reliable_mode_max = _reliable_mode_max(cosebis_nbins)
    scale_cut_full, scale_cut_fid = _scale_cuts(
        reporting, scale_cut_convention
    )

    payload = {
        "nmodes": nmodes,
        "version": version,
        "cosebis_nbins": cosebis_nbins,
        "reliable_mode_max": reliable_mode_max,
        "full_scale_cut": np.array(scale_cut_full),
        "fiducial_scale_cut": np.array(scale_cut_fid),
    }

    for name, cut in (("full", scale_cut_full), ("fiducial", scale_cut_fid)):
        print(f"Computing {name} scale cut {cut} for {version}...")
        ce_cfg, cb_cfg, cov_cfg = _compute_config_cosebis(
            xi_path, cov_path, nmodes, cut, integration
        )
        ce_h, cb_h, cov_h = _compute_harmonic_cosebis(
            pseudo_cl_path, pseudo_cov_path, nmodes, cut[0], cut[1]
        )
        payload[f"{name}_harm_En"] = ce_h
        payload[f"{name}_harm_Bn"] = cb_h
        payload[f"{name}_harm_cov"] = cov_h
        payload[f"{name}_harm_sigma_E"] = np.sqrt(
            np.clip(np.diag(cov_h[:nmodes, :nmodes]), 0, None)
        )
        payload[f"{name}_harm_sigma_B"] = np.sqrt(
            np.clip(np.diag(cov_h[nmodes:, nmodes:]), 0, None)
        )
        payload[f"{name}_config_En"] = ce_cfg
        payload[f"{name}_config_Bn"] = cb_cfg
        payload[f"{name}_config_cov"] = cov_cfg

    np.savez(out_file, **payload)
    print(f"Saved modes NPZ to {out_file}")


def _load_modes(npz_path):
    with np.load(npz_path, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


def _metrics_from_modes(modes, angular_range, out_file):
    """Emit the harmonic/config B-mode PTE + chi2 dicts for one angular range."""
    nmodes = int(modes["nmodes"])
    n_rel = int(modes["reliable_mode_max"])
    version = str(modes["version"])
    prefix = "fiducial" if angular_range == "fiducial" else "full"

    result = {}
    for path, key in (("harm", "harmonic_b_mode_ptes"), ("config", "config_b_mode_ptes")):
        bn = modes[f"{prefix}_{path}_Bn"][:n_rel]
        cov = modes[f"{prefix}_{path}_cov"][
            nmodes : nmodes + n_rel, nmodes : nmodes + n_rel
        ]
        chi2, pte, dof = compute_chi2_pte(bn, cov)
        result[key] = {
            version: {"chi2": float(chi2), "pte": float(pte), "dof": int(dof)}
        }
        print(f"  {key}: pte={pte:.4f} chi2={chi2:.2f} dof={dof}")

    result.update(
        {
            "version": version,
            "nmodes": nmodes,
            "reliable_modes": n_rel,
            "angular_range": angular_range,
            "scale_cut": [float(x) for x in modes[f"{prefix}_scale_cut"]],
            "powspace_nbins": int(modes["cosebis_nbins"]),
        }
    )

    with open(out_file, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Saved metrics to {out_file}")


def _figure_from_modes(modes, out_file):
    """Render the paper Figure 8 composite (2x2, full and fiducial cuts) from a
    modes NPZ, writing a single 300-dpi PNG."""
    nmodes = int(modes["nmodes"])
    n_rel = int(modes["reliable_mode_max"])

    def _harm(prefix):
        return {
            "En": modes[f"{prefix}_harm_En"],
            "Bn": modes[f"{prefix}_harm_Bn"],
            "sigma_E": modes[f"{prefix}_harm_sigma_E"],
            "sigma_B": modes[f"{prefix}_harm_sigma_B"],
        }

    def _cfg(prefix):
        return (
            modes[f"{prefix}_config_En"],
            modes[f"{prefix}_config_Bn"],
            modes[f"{prefix}_config_cov"],
        )

    fig = _make_combined_data_vector_figure(
        _harm("full"),
        _cfg("full"),
        _harm("fiducial"),
        _cfg("fiducial"),
        nmodes,
        tuple(float(x) for x in modes["full_scale_cut"]),
        tuple(float(x) for x in modes["fiducial_scale_cut"]),
        reliable_mode_max=n_rel,
    )
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {out_file}")


def _versions_all_for_plots():
    """Leak-corrected (non-ecut) versions plus their uncorrected counterparts."""
    leak_corr = [
        v for v in VERSIONS if "_leak_corr" in v and "_ecut" not in v
    ]
    uncorrected = [v.replace("_leak_corr", "") for v in leak_corr]
    return leak_corr + uncorrected


def _cov_integration_path(cov_dir, version, blind, min_sep, max_sep, nbins):
    """Path of the Gaussian integration-grid, masked covariance (suffix _processed.txt)."""
    base = f"covariance_{version}_{blind}_integration"
    return os.path.join(cov_dir, base, f"{base}_processed.txt")


def _angular_ranges(reporting, scale_cut_convention):
    """The full and fiducial angular ranges, by name."""
    full, fiducial = _scale_cuts(reporting, scale_cut_convention)
    return {"full": full, "fiducial": fiducial}


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(
        description="Harmonic vs config-space COSEBI cross-check (paper Figure 8)."
    )
    ap.add_argument("--fiducial-catalog", choices=FIDUCIAL_CATALOG)
    ap.add_argument("--leakage-correction", choices=["applied", "uncorrected"])
    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__)
    ap.add_argument("--integration-grid", type=INTEGRATION_GRID.__getitem__)
    ap.add_argument("--mode-count", type=MODE_COUNT.__getitem__)
    ap.add_argument(
        "--angular-range",
        choices=["full", "fiducial"],
        default="fiducial",
        help="Scale-cut regime for the per-version B-mode PTEs / primary panel",
    )
    ap.add_argument(
        "--cosmo-val-dir",
        default=None,
        help="COSMO_VAL output dir (pseudo_cl / pseudo_cl_cov FITS + xi_integration txt)",
    )
    ap.add_argument(
        "--covariance-dir",
        default=None,
        help="COSMO_INFERENCE data/covariance dir (Gaussian integration covariances)",
    )
    ap.add_argument("--blind", type=FIDUCIAL_BLIND.__getitem__)
    ap.add_argument("--patch-count", type=PATCH_COUNT.__getitem__, help="Names xi files under --cosmo-val-dir")
    ap.add_argument("--out", default=None, help="Output directory (legacy sweep mode)")
    ap.add_argument(
        "--out-file",
        default=None,
        help="Write exactly one file to this path (NPZ / JSON / PNG per mode)",
    )
    ap.add_argument(
        "--metrics-from-modes",
        default=None,
        help=(
            "Render mode: read a modes NPZ and emit the harmonic/config B-mode "
            "PTE + chi2 JSON for --angular-range to --out-file. Recomputes nothing."
        ),
    )
    ap.add_argument(
        "--figure-from-modes",
        default=None,
        help=(
            "Render mode: read a modes NPZ and render the paper Figure 8 composite "
            "(2x2, full and fiducial cuts) as a 300-dpi PNG to --out-file."
        ),
    )
    ap.add_argument(
        "--fiducial-version",
        default=None,
        help=(
            "Fiducial catalog version (paper: SP_v1.4.6.3_leak_corr). When set "
            "together with one of --fiducial-*-path below, this version's inputs "
            "are read from the given lc-reproduced path instead of the "
            "--cosmo-val-dir / --covariance-dir pattern lookup. Every other "
            "version is unaffected."
        ),
    )
    ap.add_argument(
        "--fiducial-pseudo-cl-path",
        default=None,
        help=(
            "Explicit path to the fiducial 96-bin pseudo-Cl FITS reproduced by lc "
            "(from lc's cl_bandpowers_fine; e.g. pseudo_cl_SP_v1.4.6.3_leak_corr.fits), "
            "overriding the --cosmo-val-dir pattern lookup for --fiducial-version."
        ),
    )
    ap.add_argument(
        "--fiducial-xi-path",
        default=None,
        help=(
            "Explicit path to the fiducial 1000-bin integration xi_pm text file "
            "reproduced by lc (e.g. SP_v1.4.6.3_leak_corr_xi_minsep=0.5_maxsep=300.0_"
            "nbins=1000_npatch=1.txt), overriding the --cosmo-val-dir pattern lookup "
            "for --fiducial-version."
        ),
    )
    ap.add_argument(
        "--fiducial-cov-path",
        default=None,
        help=(
            "Explicit path to the fiducial Gaussian integration covariance "
            "reproduced by lc (covariance_processed.txt), overriding the "
            "--covariance-dir pattern lookup for --fiducial-version. NOTE: there is "
            "no equivalent override for the 96-bin pseudo-Cl covariance — lc did "
            "not reproduce it (only a 32-bin cl_cov_reporting exists) — so that "
            "input always reads from the old tree, even for the fiducial version."
        ),
    )
    ap.add_argument(
        "--fiducial-pseudo-cl-cov-path",
        default=None,
        help=(
            "Explicit path to the fiducial pseudo-Cl covariance FITS matching "
            "--fiducial-pseudo-cl-path. Required by the compute pass when "
            "--cosmo-val-dir is not given."
        ),
    )
    add_scale_cut_convention_argument(ap)
    ap.add_argument("--harmonic-validation-binning", type=HARMONIC_VALIDATION_BINNING.__getitem__)
    a = ap.parse_args(argv)

    # --- Render modes: read a modes NPZ, recompute nothing, write one file. ---
    if a.metrics_from_modes or a.figure_from_modes:
        if not a.out_file:
            ap.error("--metrics-from-modes / --figure-from-modes require --out-file")
        if a.metrics_from_modes and a.figure_from_modes:
            ap.error("--metrics-from-modes and --figure-from-modes are exclusive")
        modes = _load_modes(a.metrics_from_modes or a.figure_from_modes)
        if a.metrics_from_modes:
            _metrics_from_modes(modes, a.angular_range, a.out_file)
        else:
            _figure_from_modes(modes, a.out_file)
        return

    required = ("fiducial_catalog", "leakage_correction", "harmonic_validation_binning",
                "mode_count", "integration_grid", "reporting_binning", "blind")
    missing = [f"--{r.replace('_', '-')}" for r in required if getattr(a, r) is None]
    if missing:
        ap.error(f"compute passes need {', '.join(missing)}")
    if a.cosmo_val_dir and a.patch_count is None:
        ap.error("--cosmo-val-dir needs --patch-count")
    catalogs = catalogs_from(a)
    cosebis_nbins = a.harmonic_validation_binning
    min_sep_int, max_sep_int, nbins_int = a.integration_grid
    npatch = a.patch_count

    # --- Compute pass: one version, both scale cuts, one NPZ at --out-file. ---
    if a.out_file:
        version = a.fiducial_version or catalogs.fiducial
        pseudo_cl_path = a.fiducial_pseudo_cl_path or (
            os.path.join(
                a.cosmo_val_dir,
                f"pseudo_cl_{version}_blind={a.blind}_powspace_nbins={cosebis_nbins}.fits",
            )
            if a.cosmo_val_dir
            else None
        )
        pseudo_cov_path = a.fiducial_pseudo_cl_cov_path or (
            os.path.join(
                a.cosmo_val_dir,
                f"pseudo_cl_cov_{version}_blind={a.blind}_powspace_nbins={cosebis_nbins}.fits",
            )
            if a.cosmo_val_dir
            else None
        )
        xi_path = a.fiducial_xi_path or (
            os.path.join(
                a.cosmo_val_dir,
                f"{version}_xi_minsep={min_sep_int}_maxsep={max_sep_int}"
                f"_nbins={nbins_int}_npatch={npatch}.txt",
            )
            if a.cosmo_val_dir
            else None
        )
        cov_path = a.fiducial_cov_path or (
            _cov_integration_path(
                a.covariance_dir, version, a.blind, min_sep_int, max_sep_int, nbins_int
            )
            if a.covariance_dir
            else None
        )
        missing = [
            name
            for name, value in (
                ("--fiducial-pseudo-cl-path", pseudo_cl_path),
                ("--fiducial-pseudo-cl-cov-path", pseudo_cov_path),
                ("--fiducial-xi-path", xi_path),
                ("--fiducial-cov-path", cov_path),
            )
            if value is None
        ]
        if missing:
            ap.error(
                "compute pass needs "
                + ", ".join(missing)
                + " (or --cosmo-val-dir / --covariance-dir to derive them)"
            )
        _compute_modes_npz(
            a.mode_count,
            cosebis_nbins,
            a.integration_grid,
            a.reporting_binning,
            version,
            pseudo_cl_path,
            pseudo_cov_path,
            xi_path,
            cov_path,
            a.out_file,
            a.scale_cut_convention,
        )
        return

    if not a.out:
        ap.error("one of --out-file, --out, or a --*-from-modes mode is required")

    versions = _versions_all_for_plots()
    inputs = {"specs": ["papers/bmodes/config/harmonic_config_cosebis_comparison.md"]}
    for ver in versions:
        is_fiducial = ver == a.fiducial_version

        inputs[f"pseudo_cl_{ver}"] = (
            a.fiducial_pseudo_cl_path
            if is_fiducial and a.fiducial_pseudo_cl_path
            else os.path.join(
                a.cosmo_val_dir,
                f"pseudo_cl_{ver}_blind={a.blind}_powspace_nbins={cosebis_nbins}.fits",
            )
        )
        # The 96-bin pseudo-Cl covariance has no explicit fiducial path (only the
        # 32-bin cl_cov_reporting is produced as a spec output), so every version,
        # including the fiducial, reads it from --cosmo-val-dir.
        inputs[f"pseudo_cl_cov_{ver}"] = os.path.join(
            a.cosmo_val_dir,
            f"pseudo_cl_cov_{ver}_blind={a.blind}_powspace_nbins={cosebis_nbins}.fits",
        )
        inputs[f"xi_{ver}"] = (
            a.fiducial_xi_path
            if is_fiducial and a.fiducial_xi_path
            else os.path.join(
                a.cosmo_val_dir,
                f"{ver}_xi_minsep={min_sep_int}_maxsep={max_sep_int}"
                f"_nbins={nbins_int}_npatch={npatch}.txt",
            )
        )
        inputs[f"cov_{ver}"] = (
            a.fiducial_cov_path
            if is_fiducial and a.fiducial_cov_path
            else _cov_integration_path(
                a.covariance_dir, ver, a.blind, min_sep_int, max_sep_int, nbins_int
            )
        )

    scale_cut = _angular_ranges(a.reporting_binning, a.scale_cut_convention)[a.angular_range]

    main(
        catalogs=catalogs,
        nmodes=a.mode_count,
        cosebis_nbins=cosebis_nbins,
        integration=a.integration_grid,
        reporting=a.reporting_binning,
        inputs=inputs,
        scale_cut=scale_cut,
        output_dir=a.out,
        scale_cut_convention=a.scale_cut_convention,
        paper_figure_name=f"harmonic_config_cosebis_{a.angular_range}.pdf",
        spec_path=inputs["specs"][0],
    )


if __name__ == "__main__":
    _from_cli()
