"""Raw xi_pm validation of the GLASS suite, before any E/B transform.

The campaign's other validation figures all live downstream of a pure-mode or
COSEBI filter, so a failure of the *input* -- the 20-bin xi_pm the suite ships
-- reaches them dressed as an estimator problem.  This is the undressed check:
the across-realisation mean xi_pm against the theory the mocks were generated
from ([mock_xi_input_theory]), and the across-realisation scatter against the
analytic CosmoCov sigma the paper quotes for xi_pm on the same grid.

Three rows over two columns (xi_+ | xi_-):

    1. theta xi_pm -- mock mean with the error of the mean, and the input
       theory as a line;
    2. (mean - theory) / sigma_ana -- the residual as a fraction of a SINGLE
       realisation's error bar, with sigma_emp/(sigma_ana sqrt(N)) error bars,
       so the axis reads as "fraction of an error bar" and the significance is
       the error bar rather than the point;
    3. sigma_emp / sigma_ana -- the covariance-amplitude check, with the
       +-1/sqrt(2(N-1)) sampling error of an N-realisation standard deviation.

sigma_ana is the *mock-matched* analytic sigma -- the CosmoCov run built at
the suite's own n_gal and sigma_e, rebinned to the reporting grid.  Holding
the mocks against the data covariance instead would measure the suite's 17%
galaxy excess rather than anything about the covariance, which is why
``reference_covariance`` is mockmatched throughout this campaign.  The stats
output also carries the pure shape-noise sigma the ensemble's own measured
pair counts imply, and the galaxy count those pair counts imply under a
uniform-density footprint.

Each invocation writes exactly one file:

    python mock_xi_validation.py ... --emit figure --out-file F
    python mock_xi_validation.py ... --emit stats  --out-file J
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402
import numpy as np  # noqa: E402

SRC = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC / "papers" / "bmodes" / "scripts"))
import plotting_utils  # noqa: E402
from plotting_utils import PAPER_MPLSTYLE  # noqa: E402

# Colour encodes the mode, as everywhere else in the campaign: these are the
# raw (E-dominated) correlation functions, so they take the paper's E colour,
# black, and the theory line is drawn in the same ink at lower weight.
COLOR = "black"
COLOR_THEORY = plotting_utils.CATEGORICAL_COLORS[0]

THETA_CUT = (12.0, 83.0)
SHADE = {"color": "gray", "alpha": 0.07, "lw": 0, "zorder": 0}
FIG_WIDTH = 7.1
ANNOTATION_SIZE = 12
MARK = {"ls": "none", "markersize": 4.5, "markeredgewidth": 0.8,
        "elinewidth": 0.8, "capsize": 0}

DEG2_PER_ARCMIN2 = 1.0 / 3600.0


def uniform_pair_counts(left, right, area_deg2, n_gal):
    """Pairs an annulus holds for n_gal uniform points on area_deg2.

    N_pair = (N^2 / 2) * Omega_annulus / A, the flat, uniform-density,
    no-edge prediction.  Inverting it on the measured counts is the cheapest
    probe of footprint geometry: a survey whose annuli run off the edge
    returns fewer pairs than this, and the implied N falls with theta.
    """
    annulus_deg2 = np.pi * (np.asarray(right) ** 2 - np.asarray(left) ** 2)
    annulus_deg2 = annulus_deg2 * DEG2_PER_ARCMIN2
    return 0.5 * n_gal**2 * annulus_deg2 / area_deg2


def compute(args):
    ensemble = np.load(args.reporting_ensemble, allow_pickle=False)
    theory = np.load(args.theory, allow_pickle=False)
    cov = np.load(args.covariance, allow_pickle=False)

    nbins = ensemble["xip"].shape[1]
    sigma_ana = {"xip": np.asarray(cov["sigma_xip"], dtype=float),
                 "xim": np.asarray(cov["sigma_xim"], dtype=float)}
    if sigma_ana["xip"].shape != (nbins,):
        raise SystemExit(
            f"{args.covariance} carries {sigma_ana['xip'].shape[0]} xi_+ "
            f"bins, not the ensemble's {nbins}."
        )

    theta = np.asarray(ensemble["meanr"], dtype=float).mean(axis=0)
    if not np.allclose(theta, theory["theta"], rtol=1e-6):
        raise SystemExit(
            "the theory was evaluated on a different theta grid than the "
            "ensemble reports at; it is not the matching prediction."
        )

    n_real = ensemble["xip"].shape[0]
    npairs = np.asarray(ensemble["npairs"], dtype=float).mean(axis=0)
    # Shape-noise-only sigma from the ensemble's own measured pairs.  For
    # unit weights and no intrinsic correlation, Var(xi_+) = sigma_e^4 /
    # (2 N_pair) with sigma_e the *total* dispersion (two components).
    sigma_e_total_sq = 2.0 * args.sigma_e_component**2
    sigma_shot = sigma_e_total_sq / np.sqrt(2.0 * npairs)

    implied_n = np.sqrt(
        2.0 * npairs * args.area_deg2
        / (np.pi * (ensemble["right_edges"] ** 2 - ensemble["left_edges"] ** 2)
           * DEG2_PER_ARCMIN2)
    )

    out = {
        "n_realisations": int(n_real),
        "theta_arcmin": theta.tolist(),
        "left_edges_arcmin": ensemble["left_edges"].tolist(),
        "right_edges_arcmin": ensemble["right_edges"].tolist(),
        "theta_cut_arcmin": list(THETA_CUT),
        "mean_npairs": npairs.tolist(),
        "sigma_shot_noise_xip": sigma_shot.tolist(),
        "implied_n_gal": implied_n.tolist(),
        "implied_n_gal_over_delivered": (implied_n / args.n_gal).tolist(),
        "provenance": {
            "reporting_ensemble": args.reporting_ensemble,
            "covariance": args.covariance,
            "theory": args.theory,
            "sigma_e_component": args.sigma_e_component,
            "area_deg2": args.area_deg2,
            "n_gal_delivered": args.n_gal,
            "theory_provenance": json.loads(str(theory["provenance"])),
        },
    }

    series = {}
    for key in ("xip", "xim"):
        stack = np.asarray(ensemble[key], dtype=float)
        mean = stack.mean(axis=0)
        scatter = stack.std(axis=0, ddof=1)
        err = scatter / np.sqrt(n_real)
        model = np.asarray(theory[key], dtype=float)
        model_untruncated = np.asarray(theory[key + "_untruncated"], dtype=float)
        series[key] = {
            "mean": mean, "err": err, "scatter": scatter,
            "theory": model, "theory_untruncated": model_untruncated,
            "sigma_ana": sigma_ana[key],
            "offset": (mean - model) / sigma_ana[key],
            "offset_err": scatter / sigma_ana[key] / np.sqrt(n_real),
            "z": (mean - model) / err,
            "ratio": mean / model,
            "sigma_ratio": scatter / sigma_ana[key],
        }
        out[key] = {
            "mean": mean.tolist(),
            "err_of_mean": err.tolist(),
            "theory": model.tolist(),
            "theory_untruncated": model_untruncated.tolist(),
            "mean_over_theory": (mean / model).tolist(),
            "offset_in_single_realisation_sigma": ((mean - model) / sigma_ana[key]).tolist(),
            "offset_err_in_single_realisation_sigma": (err / sigma_ana[key]).tolist(),
            "z_of_mean": ((mean - model) / err).tolist(),
            "sigma_emp": scatter.tolist(),
            "sigma_ana": sigma_ana[key].tolist(),
            "sigma_ratio": (scatter / sigma_ana[key]).tolist(),
            "sigma_emp_over_sigma_shot": (scatter / sigma_shot).tolist()
            if key == "xip" else None,
        }
    out["sigma_ratio_sampling_error"] = float(1.0 / np.sqrt(2.0 * (n_real - 1)))

    in_cut = (theta >= THETA_CUT[0]) & (theta <= THETA_CUT[1])
    out["in_cut_summary"] = {
        key: {
            "max_abs_z": float(np.max(np.abs(series[key]["z"][in_cut]))),
            "max_abs_offset_in_sigma": float(
                np.max(np.abs(series[key]["offset"][in_cut]))),
            "sigma_ratio_min": float(np.min(series[key]["sigma_ratio"][in_cut])),
            "sigma_ratio_max": float(np.max(series[key]["sigma_ratio"][in_cut])),
        }
        for key in ("xip", "xim")
    }
    return theta, series, out


def draw(theta, series, stats, out_file):
    plt.style.use(PAPER_MPLSTYLE)
    fig, axes = plt.subplots(3, 2, figsize=(FIG_WIDTH, 6.2),
                             sharex="col", sharey="row")
    ratio_err = stats["sigma_ratio_sampling_error"]

    labels = {"xip": r"$\xi_+$", "xim": r"$\xi_-$"}
    for column, key in enumerate(("xip", "xim")):
        s = series[key]
        scale = theta * 1e4
        axes[0][column].plot(
            theta, scale * s["theory_untruncated"], color=COLOR_THEORY,
            lw=1.0, ls="--", alpha=0.7, zorder=1, label=r"theory, all $\ell$")
        axes[0][column].plot(theta, scale * s["theory"], color=COLOR_THEORY,
                             lw=1.2, zorder=1,
                             label=r"theory, $\ell \leq \ell_{\rm max}^{\rm sim}$")
        axes[0][column].errorbar(
            theta, scale * s["mean"], yerr=scale * s["err"], marker="o",
            color=COLOR, ecolor=COLOR, markerfacecolor=COLOR,
            markeredgecolor=COLOR, label="mock mean", **MARK)
        axes[1][column].errorbar(
            theta, s["offset"], yerr=s["offset_err"], marker="o", color=COLOR,
            ecolor=COLOR, markerfacecolor=COLOR, markeredgecolor=COLOR, **MARK)
        axes[2][column].errorbar(
            theta, s["sigma_ratio"], yerr=ratio_err, marker="o", color=COLOR,
            ecolor=COLOR, markerfacecolor=COLOR, markeredgecolor=COLOR, **MARK)

        axes[0][column].text(
            0.04, 0.95, labels[key], transform=axes[0][column].transAxes,
            ha="left", va="top", fontsize=ANNOTATION_SIZE,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white",
                  "alpha": 0.85, "edgecolor": "none"})

    for row in range(3):
        for column in range(2):
            ax = axes[row][column]
            ax.set_xscale("log")
            ax.set_xlim(0.9, 290.0)
    for column in range(2):
        for row, level in ((1, 0.0), (2, 1.0)):
            axes[row][column].axhline(level, color="black", ls="--", lw=0.8,
                                      alpha=0.7, zorder=1)
        for row in range(3):
            ax = axes[row][column]
            low, high = ax.get_xlim()
            ax.axvspan(low, THETA_CUT[0], **SHADE)
            ax.axvspan(THETA_CUT[1], high, **SHADE)

    axes[1][0].set_ylim(-0.55, 0.55)
    axes[2][0].set_ylim(0.6, 1.4)
    axes[0][0].set_ylabel(r"$10^{4}\,\theta\,\xi\ [\mathrm{arcmin}]$")
    axes[1][0].set_ylabel(r"$(\bar{\xi} - \xi^{\rm th}) / \sigma^{\rm ana}$")
    axes[2][0].set_ylabel(r"$\sigma^{\rm emp} / \sigma^{\rm ana}$")
    for column in range(2):
        axes[2][column].set_xlabel(r"$\theta$ [arcmin]")
        ticks = [1, 10, 100]
        axes[2][column].set_xticks(ticks)
        axes[2][column].set_xticklabels([str(t) for t in ticks])
        axes[2][column].xaxis.set_minor_formatter(ticker.NullFormatter())
    axes[0][1].legend(frameon=False, loc="lower right", handletextpad=0.4,
                      borderaxespad=0.5, labelspacing=0.25)

    fig.align_ylabels(axes[:, 0])
    out = Path(out_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reporting-ensemble", required=True)
    parser.add_argument("--theory", required=True)
    parser.add_argument("--covariance", required=True,
                        help="The mock-matched reporting-grid covariance npz.")
    parser.add_argument("--sigma-e-component", type=float, default=0.26817,
                        help="The suite's per-component ellipticity dispersion.")
    parser.add_argument("--area-deg2", type=float, default=2894.0304,
                        help="Fiducial footprint area, the convention every "
                             "density in this analysis is quoted against.")
    parser.add_argument("--n-gal", type=float, default=60423575.0,
                        help="Rows per delivered realisation, for the "
                             "pair-count geometry ratio.")
    parser.add_argument("--emit", choices=("figure", "stats"), required=True)
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args()

    theta, series, stats = compute(args)
    if args.emit == "stats":
        Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_file).write_text(json.dumps(stats, indent=2) + "\n")
        print(f"Saved {args.out_file}")
    else:
        draw(theta, series, stats, args.out_file)

    for key in ("xip", "xim"):
        s = series[key]
        print(f"{key}: max |offset| = {np.max(np.abs(s['offset'])):.3f} sigma; "
              f"sigma ratio {s['sigma_ratio'].min():.3f}-"
              f"{s['sigma_ratio'].max():.3f}")


if __name__ == "__main__":
    main()
