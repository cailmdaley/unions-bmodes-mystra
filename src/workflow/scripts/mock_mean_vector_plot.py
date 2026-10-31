"""Mean B-mode data vector across GLASS mocks, per statistic.

Reads the campaign compute pass's per-realisation statistics npz -- which
carries the analytic sigma vectors alongside the empirical stack -- and plots
the across-realization mean of each B statistic in units of the analytic
(data) sigma, with sigma_emp/(sigma_ana sqrt(N)) error bars in the same
units.  The mock-to-mock scatter itself is the sigma-ratio figure's subject
and is not repeated here.
The y-axis therefore reads directly as "fraction of a data error bar":
a bias detectable at mock precision (error bars ~1/sqrt(N) ~ 0.05) can
still be negligible for the data (|mean| << 1).

Each invocation writes exactly one file:

    python mock_mean_vector_plot.py --stats S --emit figure --out-file F
    python mock_mean_vector_plot.py --stats S --emit stats  --out-file J
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import yaml

# The two script trees are separate sys.path roots; the fiducial window is
# resolved by the paper-facing helper rather than restated here.
sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from scale_cuts import (  # noqa: E402
    add_scale_cut_convention_argument,
    resolve_fiducial_bin_window,
    resolve_fiducial_scale_cuts,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

try:
    from plotting_utils import (  # type: ignore
        CATEGORICAL_COLORS,
        PAPER_MPLSTYLE,
        SCATTER_ERRORBAR,
    )

    plt.style.use(PAPER_MPLSTYLE)
except Exception:  # pragma: no cover - style is cosmetic
    PAPER_MPLSTYLE = None
    CATEGORICAL_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
    SCATTER_ERRORBAR = {"ls": "none", "ms": 4, "elinewidth": 0.8, "capsize": 0}


def _hartlap(n_real: int, p: int) -> float:
    return (n_real - p - 2) / (n_real - 1)


def _mean_chi2_pte(vectors: np.ndarray) -> tuple[float, float, int]:
    """Chi2/PTE of the mean against the empirical covariance of the mean."""
    from scipy import stats

    n_real, p = vectors.shape
    mean = vectors.mean(axis=0)
    cov_mean = np.cov(vectors, rowvar=False, ddof=1) / n_real
    inv = np.linalg.inv(cov_mean) * _hartlap(n_real, p)
    chi2 = float(mean @ inv @ mean)
    return chi2, float(stats.chi2.sf(chi2, p)), p


def _scale_averaged_stats(
    vectors: np.ndarray, sigma_ana: np.ndarray, bins: np.ndarray
) -> dict[str, object]:
    """Summarize per-mock scale averages and the largest mean bias."""
    normalized = vectors / sigma_ana
    results: dict[str, object] = {}
    for name, values in (("raw", vectors), ("normalized", normalized)):
        per_mock = values.mean(axis=1)
        mean = float(per_mock.mean())
        std = float(per_mock.std(ddof=1))
        stderr = std / np.sqrt(vectors.shape[0])
        results[name] = {
            "mean": mean,
            "std": std,
            "stderr": float(stderr),
            "z": float(mean / stderr),
        }

    mean_normalized = normalized.mean(axis=0)
    max_index = int(np.argmax(np.abs(mean_normalized)))
    results["max_abs_mean_normalized"] = float(abs(mean_normalized[max_index]))
    results["max_abs_mean_normalized_bin"] = bins[max_index].item()
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--emit", choices=("figure", "stats"), required=True)
    parser.add_argument("--out-file", required=True)
    parser.add_argument("--label", default=None)
    parser.add_argument("--config", required=True, help="bmodes config.yaml")
    add_scale_cut_convention_argument(parser)
    args = parser.parse_args()

    with open(args.config) as stream:
        config = yaml.safe_load(stream)
    # The pure-mode statistics are quoted on two ranges, and the file says
    # which is which: the full 20-bin reporting grid, and the fiducial window
    # the null test is actually claimed on. They differ materially — the
    # mocks carry a small mean B residual at the smallest scales, which the
    # cut excludes — so an unlabelled single number invites the two to be
    # read against each other.
    window = resolve_fiducial_bin_window(config, args.scale_cut_convention)
    window_arcmin = resolve_fiducial_scale_cuts(config, args.scale_cut_convention)

    data = np.load(args.stats, allow_pickle=False)
    provenance = json.loads(str(data["provenance"]))
    label = args.label or provenance.get("covariance_mode", {}).get(
        "label", Path(args.stats).stem
    )
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_real = data["Bn"].shape[0]

    panels = [
        (
            "Bn",
            np.arange(1, 21),
            np.asarray(data["sigma_cosebis"]),
            "COSEBI mode $n$",
            r"$\langle B_n\rangle/\sigma_n^\mathrm{ana}$",
            None,
            [("$B_n$", "Bn", "o")],
        ),
        (
            "pure",
            np.asarray(data["theta"]),
            None,
            r"$\theta$ [arcmin]",
            r"$\langle\xi^B\rangle/\sigma^\mathrm{ana}$",
            "log",
            [
                (r"$\xi_+^B$", "xip_B_full", "o"),
                (r"$\xi_-^B$", "xim_B_full", "s"),
            ],
        ),
        (
            "cl",
            np.asarray(data["ell_cut"]),
            np.asarray(data["sigma_cl_cut"]),
            r"$\ell$",
            r"$\langle C_\ell^{BB}\rangle/\sigma_\ell^\mathrm{ana}$",
            None,
            [(r"$C_\ell^{BB}$", "cl_BB", "o")],
        ),
    ]
    pure_sigma = {
        "xip_B_full": np.asarray(data["sigma_xip_full"]),
        "xim_B_full": np.asarray(data["sigma_xim_full"]),
    }

    fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.8))
    all_stats: dict[str, dict[str, object]] = {}
    for ax, (name, x, sigma_ana, xlabel, ylabel, xscale, series) in zip(axes, panels):
        annotations = []
        for colour, (series_label, key, marker) in zip(CATEGORICAL_COLORS, series):
            vectors = np.asarray(data[key])
            sig_ana = pure_sigma[key] if sigma_ana is None else sigma_ana
            mean = vectors.mean(axis=0) / sig_ana
            scatter = vectors.std(axis=0, ddof=1) / sig_ana
            err = scatter / np.sqrt(n_real)
            ax.errorbar(
                x,
                mean,
                yerr=err,
                marker=marker,
                color=colour,
                ecolor=colour,
                label=series_label,
                **SCATTER_ERRORBAR,
            )
            chi2, pte, dof = _mean_chi2_pte(vectors)
            stat_name = {
                "Bn": "B_n",
                "xip_B_full": "xi_p_B",
                "xim_B_full": "xi_m_B",
                "cl_BB": "C_ell_BB",
            }[key]
            entry: dict[str, object] = {
                "n_bins": int(vectors.shape[1]),
                "bins": [value.item() for value in x],
                "scale_averaged": _scale_averaged_stats(vectors, sig_ana, x),
                "full_range": {
                    "chi2": chi2,
                    "pte": pte,
                    "dof": dof,
                    "n_bins": int(vectors.shape[1]),
                    "range": "all reporting bins",
                },
            }
            annotations.append(
                rf"{series_label}: $\chi^2$={chi2:.1f}/{dof}, PTE={pte:.3g}"
            )
            if name == "pure":
                lo, hi = window
                cut_chi2, cut_pte, cut_dof = _mean_chi2_pte(vectors[:, lo : hi + 1])
                entry["fiducial_window"] = {
                    "chi2": cut_chi2,
                    "pte": cut_pte,
                    "dof": cut_dof,
                    "n_bins": cut_dof,
                    "bin_indices": [lo, hi],
                    "theta_range_arcmin": list(window_arcmin),
                    "scale_cut_convention": args.scale_cut_convention,
                    "range": "fiducial scale cut",
                }
                annotations.append(
                    rf"    cut: $\chi^2$={cut_chi2:.1f}/{cut_dof}, "
                    rf"PTE={cut_pte:.3g}"
                )
            all_stats[stat_name] = entry
        ax.axhline(0.0, color="black", ls="--", lw=0.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xscale:
            ax.set_xscale(xscale)
        legend_handles = [
            Line2D(
                [],
                [],
                color="0.35",
                lw=1.0,
                label=r"$\sigma_\mathrm{emp}/(\sigma_\mathrm{ana}\sqrt{N})$",
            ),
        ]
        handles, _ = ax.get_legend_handles_labels()
        legend_handles = handles + legend_handles
        ax.legend(handles=legend_handles, fontsize=5.5, loc="lower right")
        ax.text(
            0.03,
            0.97,
            "\n".join(annotations),
            transform=ax.transAxes,
            fontsize=5,
            va="top",
        )
        ax.grid(alpha=0.2)
    axes[0].set_xticks([1, 5, 10, 15, 20])
    # One y-range across the three panels, so a mean that is large in one
    # statistic reads as large next to the others rather than being rescaled
    # away by its own panel's autoscale.
    limits = [axis.get_ylim() for axis in axes]
    shared = (min(low for low, _ in limits), max(high for _, high in limits))
    for axis in axes:
        axis.set_ylim(shared)
    fig.suptitle(
        rf"GLASS mock campaign ({label}): mean B vector in units of analytic $\sigma$"
        rf" (error bars: $\sigma_\mathrm{{emp}}/("
        rf"\sigma_\mathrm{{ana}}\sqrt{{N}})$, N={n_real})",
        fontsize=8,
    )
    fig.tight_layout()
    if args.emit == "figure":
        fig.savefig(out_path, dpi=250, bbox_inches="tight")
        print(f"Saved {out_path}")
    else:
        with open(out_path, "w") as handle:
            json.dump(
                {
                    "n_realizations": n_real,
                    "fiducial_window": {
                        "bin_indices": list(window),
                        "theta_range_arcmin": list(window_arcmin),
                        "scale_cut_convention": args.scale_cut_convention,
                    },
                    "statistics": all_stats,
                },
                handle,
                indent=2,
            )
            handle.write("\n")
        print(f"Saved {out_path}")
    plt.close(fig)

    for name, values in all_stats.items():
        full = values["full_range"]
        cut = values.get("fiducial_window")
        scale = values["scale_averaged"]
        raw = scale["raw"]
        normalized = scale["normalized"]
        print(
            f"mean-vector {name}: "
            f"raw mean={raw['mean']:.6g}, std={raw['std']:.6g}, "
            f"stderr={raw['stderr']:.6g}, z={raw['z']:.6g}; "
            f"normalized mean={normalized['mean']:.6g}, "
            f"std={normalized['std']:.6g}, stderr={normalized['stderr']:.6g}, "
            f"z={normalized['z']:.6g}; "
            f"max |mean|/sigma_ana={scale['max_abs_mean_normalized']:.6g} "
            f"(bin {scale['max_abs_mean_normalized_bin']}); "
            f"full-range chi2={full['chi2']:.1f}/{full['dof']} "
            f"PTE={full['pte']:.4g}"
            + (
                f"; fiducial-window chi2={cut['chi2']:.1f}/{cut['dof']} "
                f"PTE={cut['pte']:.4g}"
                if cut
                else ""
            )
        )


if __name__ == "__main__":
    main()
