"""
BB Covariance Blind Independence Claim

Tests whether BB covariances are blind-independent (as expected for null signals)
while EE covariances vary across blinds (due to sample variance from cosmological signal).

Compares:
- Pure E/B: cov(xi+^B), cov(xi-^B) vs cov(xi+^E), cov(xi-^E)
- COSEBIS: cov(B_n) vs cov(E_n)
- Harmonic: cov(C_ell^BB) vs cov(C_ell^EE)
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import treecorr
from astropy.io import fits
from plotting_utils import PAPER_MPLSTYLE

from sp_validation.b_modes import calculate_cosebis
from decision_values import INTEGRATION_GRID, MC_SAMPLES, MODE_COUNT, REPORTING_BINNING

plt.style.use(PAPER_MPLSTYLE)


BLINDS = ["A", "B", "C"]


def load_pure_eb_diagonals(path):
    """Load pure E/B covariance and extract block diagonals.

    Covariance is 120x120 with 6 blocks of 20 bins:
    [E+, E-, B+, B-, amb+, amb-]
    """
    data = np.load(path)
    cov = data["cov_pure_eb"]
    theta = data["theta"]
    nbins = len(theta)

    def get_block_diag(block_idx):
        sl = slice(block_idx * nbins, (block_idx + 1) * nbins)
        return np.diag(cov[sl, sl])

    return {
        "theta": theta,
        "xip_E": get_block_diag(0),
        "xim_E": get_block_diag(1),
        "xip_B": get_block_diag(2),
        "xim_B": get_block_diag(3),
    }


def load_harmonic_diagonals(path):
    """Load harmonic covariance and extract EE and BB diagonals."""
    with fits.open(path) as hdu:
        diag_EE = np.diag(hdu["COVAR_EE_EE"].data)
        diag_BB = np.diag(hdu["COVAR_BB_BB"].data)
    return {
        "EE": diag_EE,
        "BB": diag_BB,
    }


def load_cosebis_diagonals(
    xi_integration_path,
    cov_integration_path,
    nmodes,
    theta_min,
    theta_max,
    min_sep_int,
    max_sep_int,
    nbins_int,
):
    """Compute COSEBIS covariance diagonals from config-space covariance.

    Uses calculate_cosebis to transform config-space covariance to COSEBIS space.
    Returns E_n and B_n covariance diagonals.
    """
    # Load fine-binned 2PCF (need the binning info for COSEBIS calculation)
    gg = treecorr.GGCorrelation(
        min_sep=min_sep_int, max_sep=max_sep_int, nbins=nbins_int, sep_units="arcmin"
    )
    gg.read(xi_integration_path)

    # Compute COSEBIS with this blind's covariance
    results = calculate_cosebis(
        gg,
        nmodes=nmodes,
        scale_cuts=[(theta_min, theta_max)],
        cov_path=cov_integration_path,
    )

    # Extract covariance for the fiducial scale cut
    result = results[(theta_min, theta_max)]
    cov = result["cov"]

    # E is [:nmodes, :nmodes], B is [nmodes:, nmodes:]
    diag_E = np.diag(cov[:nmodes, :nmodes])
    diag_B = np.diag(cov[nmodes:, nmodes:])

    return {
        "E": diag_E,
        "B": diag_B,
        "nmodes": nmodes,
    }


def compute_ratios(diag_ref, diag_test):
    """Compute ratio and deviation statistics."""
    ratio = diag_test / diag_ref
    dev = np.abs(ratio - 1.0)
    return {
        "ratio": ratio,
        "max_dev": float(np.max(dev)),
        "mean_dev": float(np.mean(dev)),
        "min_ratio": float(np.min(ratio)),
        "max_ratio": float(np.max(ratio)),
    }


def make_figure(
    theta,
    ell_eff,
    pure_eb_results,
    harmonic_results,
    cosebis_results,
    output_path,
    n_samples=2000,
):
    """Four-panel figure comparing BB vs EE stability across blinds.

    Layout: 2x2
    - Top row: Pure E/B (xi+, xi-)
    - Bottom row: COSEBIS (B_n vs E_n), Harmonic (C_ell)

    Color encoding: E-mode vs B-mode
    Marker encoding: square = B/A, triangle = C/A (same for both E and B)
    """
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))

    # Colors: E-mode vs B-mode
    color_E = "#E69F00"  # orange for E
    color_B = "#0072B2"  # blue for B

    # Expected 1σ error on ratio of two MC-estimated quantities
    # σ(ratio) ≈ √(2/N) for ratio ≈ 1
    ratio_err = np.sqrt(2.0 / n_samples)

    # Helper for ratio panels
    def setup_ratio_panel(ax, ylabel_x, title, show_mc_band=False):
        ax.axhline(1.0, color="gray", ls="-", lw=0.8, zorder=0)
        if show_mc_band:
            # Show expected 1σ MC scatter band only
            ax.axhspan(
                1 - ratio_err,
                1 + ratio_err,
                color="gray",
                alpha=0.25,
                label=rf"$\pm\sqrt{{2/N}}$ ($N={n_samples}$)",
            )
        ax.set_xscale("log")
        ax.set_xlabel(ylabel_x)
        ax.set_title(title)

    # --- Panel 1: xi+^B vs xi+^E ---
    ax = axes[0, 0]
    setup_ratio_panel(
        ax,
        r"$\theta$ [arcmin]",
        r"$\xi_+$: covariance ratio across blinds",
        show_mc_band=True,
    )

    # B-mode (blue): square for B/A, triangle for C/A
    ax.errorbar(
        theta,
        pure_eb_results["B"]["xip_B"]["ratio"],
        yerr=ratio_err,
        fmt="s",
        color=color_B,
        label=r"B-mode B/A",
        markersize=5,
        alpha=0.8,
        capsize=0,
    )
    ax.errorbar(
        theta * 1.03,
        pure_eb_results["C"]["xip_B"]["ratio"],
        yerr=ratio_err,
        fmt="^",
        color=color_B,
        label=r"B-mode C/A",
        markersize=5,
        alpha=0.8,
        capsize=0,
    )

    # E-mode (orange): square for B/A, triangle for C/A
    ax.plot(
        theta,
        pure_eb_results["B"]["xip_E"]["ratio"],
        "s",
        color=color_E,
        label=r"E-mode B/A",
        markersize=4,
        alpha=0.6,
    )
    ax.plot(
        theta * 1.03,
        pure_eb_results["C"]["xip_E"]["ratio"],
        "^",
        color=color_E,
        label=r"E-mode C/A",
        markersize=4,
        alpha=0.6,
    )

    ax.set_ylabel("Diagonal ratio")
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.set_xlim(1, 300)
    ax.set_ylim(0.85, 1.15)

    # --- Panel 2: xi-^B vs xi-^E ---
    ax = axes[0, 1]
    setup_ratio_panel(
        ax,
        r"$\theta$ [arcmin]",
        r"$\xi_-$: covariance ratio across blinds",
        show_mc_band=True,
    )

    ax.errorbar(
        theta,
        pure_eb_results["B"]["xim_B"]["ratio"],
        yerr=ratio_err,
        fmt="s",
        color=color_B,
        label=r"B-mode B/A",
        markersize=5,
        alpha=0.8,
        capsize=0,
    )
    ax.errorbar(
        theta * 1.03,
        pure_eb_results["C"]["xim_B"]["ratio"],
        yerr=ratio_err,
        fmt="^",
        color=color_B,
        label=r"B-mode C/A",
        markersize=5,
        alpha=0.8,
        capsize=0,
    )

    ax.plot(
        theta,
        pure_eb_results["B"]["xim_E"]["ratio"],
        "s",
        color=color_E,
        label=r"E-mode B/A",
        markersize=4,
        alpha=0.6,
    )
    ax.plot(
        theta * 1.03,
        pure_eb_results["C"]["xim_E"]["ratio"],
        "^",
        color=color_E,
        label=r"E-mode C/A",
        markersize=4,
        alpha=0.6,
    )

    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.set_xlim(1, 300)
    ax.set_ylim(0.85, 1.15)

    # --- Panel 3: COSEBIS B_n vs E_n ---
    ax = axes[1, 0]
    nmodes = cosebis_results["B"]["nmodes"]
    n_arr = np.arange(1, nmodes + 1)
    ax.axhline(1.0, color="gray", ls="-", lw=0.8, zorder=0)
    ax.set_xlabel(r"Mode $n$")
    ax.set_title(r"COSEBIS: covariance ratio across blinds")

    # B-mode (blue): square for B/A, triangle for C/A
    ax.plot(
        n_arr,
        cosebis_results["B"]["B"]["ratio"],
        "s",
        color=color_B,
        label=r"B-mode B/A",
        markersize=5,
        alpha=0.8,
    )
    ax.plot(
        n_arr + 0.15,
        cosebis_results["C"]["B"]["ratio"],
        "^",
        color=color_B,
        label=r"B-mode C/A",
        markersize=5,
        alpha=0.8,
    )

    # E-mode (orange): square for B/A, triangle for C/A
    ax.plot(
        n_arr,
        cosebis_results["B"]["E"]["ratio"],
        "s",
        color=color_E,
        label=r"E-mode B/A",
        markersize=4,
        alpha=0.6,
    )
    ax.plot(
        n_arr + 0.15,
        cosebis_results["C"]["E"]["ratio"],
        "^",
        color=color_E,
        label=r"E-mode C/A",
        markersize=4,
        alpha=0.6,
    )

    ax.set_ylabel("Diagonal ratio")
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.set_ylim(0.85, 1.15)

    # --- Panel 4: C_ell^BB vs C_ell^EE ---
    ax = axes[1, 1]
    setup_ratio_panel(ax, r"$\ell$", r"$C_\ell$: covariance ratio across blinds")

    # B-mode (blue): square for B/A, triangle for C/A
    ax.plot(
        ell_eff,
        harmonic_results["B"]["BB"]["ratio"],
        "s",
        color=color_B,
        label=r"BB B/A",
        markersize=5,
        alpha=0.8,
    )
    ax.plot(
        ell_eff * 1.03,
        harmonic_results["C"]["BB"]["ratio"],
        "^",
        color=color_B,
        label=r"BB C/A",
        markersize=5,
        alpha=0.8,
    )

    # E-mode (orange): square for B/A, triangle for C/A
    ax.plot(
        ell_eff,
        harmonic_results["B"]["EE"]["ratio"],
        "s",
        color=color_E,
        label=r"EE B/A",
        markersize=4,
        alpha=0.6,
    )
    ax.plot(
        ell_eff * 1.03,
        harmonic_results["C"]["EE"]["ratio"],
        "^",
        color=color_E,
        label=r"EE C/A",
        markersize=4,
        alpha=0.6,
    )

    ax.set_ylabel("Diagonal ratio")
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.set_ylim(0.85, 1.15)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(
    integration,
    n_samples,
    pure_eb_paths,
    harmonic_paths,
    xi_integration_path,
    cov_integration_paths,
    pseudo_cl_path,
    nmodes,
    theta_min,
    theta_max,
    figure_path,
):
    """Run the BB-covariance blind-independence cross-check.

    ``pure_eb_paths`` / ``harmonic_paths`` / ``cov_integration_paths`` are dicts
    keyed by blind (A/B/C). All paths are absolute; the per-blind mock-version
    covariances are read directly, so the check is self-contained.
    """
    # Load pure E/B covariances for all blinds
    pure_eb_data = {}
    for blind in BLINDS:
        pure_eb_data[blind] = load_pure_eb_diagonals(pure_eb_paths[blind])

    theta = pure_eb_data["A"]["theta"]

    # Load harmonic covariances for all blinds
    harmonic_data = {}
    for blind in BLINDS:
        harmonic_data[blind] = load_harmonic_diagonals(harmonic_paths[blind])

    min_sep_int, max_sep_int, nbins_int = integration

    cosebis_data = {}
    for blind in BLINDS:
        cosebis_data[blind] = load_cosebis_diagonals(
            xi_integration_path,
            cov_integration_paths[blind],
            nmodes,
            theta_min,
            theta_max,
            min_sep_int,
            max_sep_int,
            nbins_int,
        )

    # Read ell bin centers from pseudo-Cl data file
    with fits.open(pseudo_cl_path) as hdu:
        ell_eff = hdu["PSEUDO_CELL"].data["ELL"]

    # Compute ratios relative to blind A
    pure_eb_results = {}
    for blind in ["B", "C"]:
        pure_eb_results[blind] = {}
        for mode in ["xip_E", "xim_E", "xip_B", "xim_B"]:
            pure_eb_results[blind][mode] = compute_ratios(
                pure_eb_data["A"][mode], pure_eb_data[blind][mode]
            )

    harmonic_results = {}
    for blind in ["B", "C"]:
        harmonic_results[blind] = {}
        for mode in ["EE", "BB"]:
            harmonic_results[blind][mode] = compute_ratios(
                harmonic_data["A"][mode], harmonic_data[blind][mode]
            )

    cosebis_results = {}
    for blind in ["B", "C"]:
        cosebis_results[blind] = {
            "nmodes": nmodes,
        }
        for mode in ["E", "B"]:
            cosebis_results[blind][mode] = compute_ratios(
                cosebis_data["A"][mode], cosebis_data[blind][mode]
            )

    # Generate figure
    make_figure(
        theta,
        ell_eff,
        pure_eb_results,
        harmonic_results,
        cosebis_results,
        figure_path,
        n_samples=n_samples,
    )

    # Compute summary statistics
    bb_max_devs = [
        pure_eb_results["B"]["xip_B"]["max_dev"],
        pure_eb_results["C"]["xip_B"]["max_dev"],
        pure_eb_results["B"]["xim_B"]["max_dev"],
        pure_eb_results["C"]["xim_B"]["max_dev"],
        cosebis_results["B"]["B"]["max_dev"],
        cosebis_results["C"]["B"]["max_dev"],
        harmonic_results["B"]["BB"]["max_dev"],
        harmonic_results["C"]["BB"]["max_dev"],
    ]
    ee_max_devs = [
        pure_eb_results["B"]["xip_E"]["max_dev"],
        pure_eb_results["C"]["xip_E"]["max_dev"],
        pure_eb_results["B"]["xim_E"]["max_dev"],
        pure_eb_results["C"]["xim_E"]["max_dev"],
        cosebis_results["B"]["E"]["max_dev"],
        cosebis_results["C"]["E"]["max_dev"],
        harmonic_results["B"]["EE"]["max_dev"],
        harmonic_results["C"]["EE"]["max_dev"],
    ]

    # Pass criteria differ between the analytic and MC methods:
    # COSEBIS uses analytic propagation (T @ Cov @ T.T) → blind-independent
    # Pure E/B and Harmonic use MC sampling → inherit sampling noise
    cosebis_bb_max = max(
        cosebis_results["B"]["B"]["max_dev"],
        cosebis_results["C"]["B"]["max_dev"],
    )
    cosebis_ee_max = max(
        cosebis_results["B"]["E"]["max_dev"],
        cosebis_results["C"]["E"]["max_dev"],
    )

    # COSEBIS B_n should be blind-independent (<0.1%)
    cosebis_bb_blind_independent = cosebis_bb_max < 0.001
    # E-modes should vary ~10% due to sample variance
    ee_varies_as_expected = 0.05 < cosebis_ee_max < 0.15

    # Summary over all three statistics
    bb_max = max(bb_max_devs)
    ee_max = max(ee_max_devs)
    bb_closer_to_unity = bb_max < ee_max
    bb_within_2pct = bb_max < 0.02

    # Print summary
    print("\nBB Covariance Blind Independence Summary:")
    print(f"  COSEBIS B_n max deviation: {cosebis_bb_max * 100:.6f}%")
    print(f"  COSEBIS E_n max deviation: {cosebis_ee_max * 100:.2f}%")
    print(f"  COSEBIS B_n blind-independent (<0.1%): {cosebis_bb_blind_independent}")
    print(f"  E-modes vary as expected (5-15%): {ee_varies_as_expected}")
    print("\n  Pure E/B + Harmonic (MC methods):")
    print(f"    BB max deviation: {bb_max * 100:.2f}%")
    print(f"    EE max deviation: {ee_max * 100:.2f}%")
    print("    Note: MC sampling noise causes BB to vary similarly to EE")


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(
        description="BB-covariance blind-independence cross-check (mock version)."
    )
    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)
    ap.add_argument("--integration-grid", type=INTEGRATION_GRID.__getitem__, required=True)
    ap.add_argument("--mode-count", type=MODE_COUNT.__getitem__, required=True)
    ap.add_argument("--mc-samples", type=MC_SAMPLES.__getitem__, required=True)
    versions = {
        "mock_version_v146": "SP_v1.4.6_leak_corr",
        "size_cut_v1463": "SP_v1.4.6.3_leak_corr",
    }
    ap.add_argument("--version", required=True, choices=versions)
    ap.add_argument("--pure-eb-template", required=True)
    ap.add_argument("--harmonic-template", required=True)
    ap.add_argument("--cov-integration-template", required=True)
    ap.add_argument("--xi-template", required=True)
    ap.add_argument("--pseudo-cl-template", required=True)
    ap.add_argument("--out-file", required=True, help="Output PNG path (lc {output})")
    a = ap.parse_args(argv)

    version = versions[a.version]
    pure_eb_paths = {b: a.pure_eb_template.format(version=version, blind=b) for b in BLINDS}
    harmonic_paths = {b: a.harmonic_template.format(version=version, blind=b) for b in BLINDS}
    cov_integration_paths = {
        b: a.cov_integration_template.format(version=version, blind=b) for b in BLINDS
    }
    xi_integration_path = a.xi_template.format(version=version)
    pseudo_cl_path = a.pseudo_cl_template.format(version=version)

    main(
        integration=a.integration_grid,
        n_samples=a.mc_samples,
        pure_eb_paths=pure_eb_paths,
        harmonic_paths=harmonic_paths,
        xi_integration_path=xi_integration_path,
        cov_integration_paths=cov_integration_paths,
        pseudo_cl_path=pseudo_cl_path,
        nmodes=a.mode_count,
        theta_min=a.reporting_binning.min_sep,
        theta_max=a.reporting_binning.max_sep,
        figure_path=a.out_file,
    )


if __name__ == "__main__":
    _from_cli()
