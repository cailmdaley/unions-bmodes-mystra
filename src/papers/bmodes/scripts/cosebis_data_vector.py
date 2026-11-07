"""COSEBI E_n and B_n, and the paper's B_n / sigma_n figure.

Two single-file modes:
- compute: --xi-integration --out-file <npz>                             (E_n/B_n at both cuts)
- render:  --plot-from-modes <npz> --covariance <npz> --out-file <png>   (paper figure, 300 dpi)
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import treecorr
from matplotlib.lines import Line2D
from plotting_utils import (
    COLOR_COSEBIS,
    COLOR_COSEBIS_ALL,
    DATA_VECTOR_MARK,
    EXCLUDED_SHADE,
    FIG_WIDTH_SINGLE,
    PAPER_MPLSTYLE,
)
from print_width import savefig_at_width
from scale_cuts import add_scale_cut_convention_argument, resolve_fiducial_scale_cuts

from cosebis_transform import cosebis_modes
from decision_values import INTEGRATION_GRID, MODE_COUNT, REPORTING_BINNING

plt.style.use(PAPER_MPLSTYLE)


def _create_single_panel_bmode_figure(datasets, nmodes, scale_cuts, title=None):
    """Create a single-panel B-mode COSEBIS figure.

    Plots B_n / sigma_n (dimensionless, in units of standard deviation).
    The COSEBI family hue (purple): filled circles in the saturated hue for
    the fiducial cuts, open diamonds in the light tint for the full range,
    as in the mock-validation figure.
    Modes beyond the n <= 6 subset are shaded, the way scales outside the
    fiducial cuts are shaded elsewhere.

    Args:
        title: Optional title for the figure (None for paper figure).
    """
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_SINGLE, FIG_WIDTH_SINGLE * 0.47))

    modes = np.arange(1, nmodes + 1)
    # scale key -> (x offset, marker, filled, colour, name)
    series = {
        "fiducial": (-0.15, "o", True, COLOR_COSEBIS, "fiducial cuts"),
        "full": (0.15, "D", False, COLOR_COSEBIS_ALL, "all scales"),
    }

    handles = []
    for scale_key, data in datasets.items():
        offset, marker, filled, color, name = series[scale_key]
        lo, hi = scale_cuts[scale_key]
        label = rf"{name} (${lo:.0f}$--${hi:.0f}'$)"
        face = color if filled else "none"
        ax.errorbar(
            modes + offset,
            data["Bn_normalized"],
            yerr=np.ones_like(data["Bn_normalized"]),
            marker=marker,
            color=color,
            markerfacecolor=face,
            markeredgecolor=color,
            **DATA_VECTOR_MARK,
        )
        # Marker-only legend entries: the bars are all unit length by
        # construction and would only add clutter.
        handles.append(
            Line2D(
                [], [], ls="none", marker=marker, color=color,
                markerfacecolor=face, markeredgecolor=color,
                markeredgewidth=DATA_VECTOR_MARK["markeredgewidth"],
                markersize=DATA_VECTOR_MARK["markersize"], label=label,
            )
        )

    ax.axhline(0, color="black", ls="--", linewidth=0.8, alpha=0.6, zorder=1)
    ax.axvspan(6.5, nmodes + 0.5, **EXCLUDED_SHADE)
    ax.set_ylabel(r"$B_n / \sigma_n$")
    ax.set_xlabel("COSEBI mode $n$")
    ax.set_xlim(0.5, nmodes + 0.5)
    ax.set_xticks(np.arange(1, nmodes + 1))
    ax.set_xticklabels(
        [str(i) for i in range(1, nmodes + 1)],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    ax.tick_params(axis="both", width=0.5, length=3)

    if title:
        ax.set_title(f"COSEBI B-modes ({title})")

    # Unit error bars span roughly -2.8 to 5.1; leave room for the legend.
    ax.set_ylim(-3.5, 5.5)
    ax.set_yticks(np.arange(-3, 6, 2))

    ax.legend(
        handles=handles,
        loc="upper right",
        ncol=2,
        frameon=False,
        handletextpad=0.3,
        columnspacing=0.8,
        borderaxespad=0.3,
    )

    plt.tight_layout()
    return fig


def _plot_from_modes(modes_npz, covariance_npz, out_file):
    """Render the paper data-vector figure from the modes and their covariance."""
    d = np.load(modes_npz, allow_pickle=True)
    cov = np.load(covariance_npz, allow_pickle=True)
    nmodes = int(d["nmodes"])
    scale_cuts = {}
    datasets = {}
    for scale_key in ("fiducial", "full"):
        scale_cuts[scale_key] = tuple(float(v) for v in d[f"{scale_key}_scale_cut"])
        datasets[scale_key] = {
            "Bn_normalized": d[f"{scale_key}_Bn"] / cov[f"{scale_key}_sigma_B"]
        }
    fig = _create_single_panel_bmode_figure(datasets, nmodes, scale_cuts, title=None)
    savefig_at_width(fig, out_file, FIG_WIDTH_SINGLE, dpi=300)
    plt.close(fig)
    print(f"Saved figure to {out_file}")


def main(version, nmodes, integration, reporting, xi_integration, out_file, scale_cut_convention):

    scale_cuts = {
        "fiducial": resolve_fiducial_scale_cuts(reporting, scale_cut_convention),
        "full": (reporting.min_sep, reporting.max_sep),
    }

    min_sep_int, max_sep_int, nbins_int = integration

    gg = treecorr.GGCorrelation(
        min_sep=min_sep_int,
        max_sep=max_sep_int,
        nbins=nbins_int,
        sep_units="arcmin",
    )
    gg.read(xi_integration)

    npz_payload = {"nmodes": nmodes, "version": version}
    for scale_key, scale_cut in scale_cuts.items():
        En, Bn = cosebis_modes(gg, nmodes, scale_cut)
        npz_payload[f"{scale_key}_En"] = En
        npz_payload[f"{scale_key}_Bn"] = Bn
        npz_payload[f"{scale_key}_scale_cut"] = np.array(scale_cut)

    with open(out_file, "wb") as f:
        np.savez(f, **npz_payload)
    print(f"Saved COSEBI modes data to {out_file}")


def _from_cli(argv=None):

    ap = argparse.ArgumentParser(
        description="COSEBI E_n/B_n data vector (NPZ) and paper figure."
    )
    ap.add_argument("--version", help="Catalog version the modes are measured on")
    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__)
    ap.add_argument("--integration-grid", type=INTEGRATION_GRID.__getitem__)
    ap.add_argument("--mode-count", type=MODE_COUNT.__getitem__)
    ap.add_argument(
        "--xi-integration",
        help="Fine-binned (1000-bin) TreeCorr xi_pm .txt dump",
    )
    ap.add_argument(
        "--plot-from-modes",
        help="Plot-only mode: render the paper figure from a modes NPZ",
    )
    ap.add_argument(
        "--covariance",
        help="Plot-only mode: the COSEBI covariance NPZ that sets sigma_n",
    )
    ap.add_argument(
        "--out-file",
        required=True,
        help="Output file (lc {output}): the modes NPZ, or the PNG in --plot-from-modes mode",
    )
    add_scale_cut_convention_argument(ap)
    a = ap.parse_args(argv)
    if a.plot_from_modes:
        _plot_from_modes(a.plot_from_modes, a.covariance, a.out_file)
        return
    required = ("xi_integration", "version", "reporting_binning", "integration_grid", "mode_count")
    missing = [f"--{r.replace('_', '-')}" for r in required if getattr(a, r) is None]
    if missing:
        ap.error(f"compute mode needs {', '.join(missing)}")
    main(a.version, a.mode_count, a.integration_grid, a.reporting_binning,
         a.xi_integration, a.out_file, a.scale_cut_convention)


if __name__ == "__main__":
    _from_cli()
