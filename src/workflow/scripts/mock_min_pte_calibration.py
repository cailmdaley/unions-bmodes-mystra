#!/usr/bin/env python3
"""Calibrate a global minimum-PTE criterion from the GLASS mocks."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np

# The two script trees are separate sys.path roots; the fiducial window is
# resolved by the paper-facing helper rather than restated here, so this
# calibration and the paper's PTE matrices cut on the same bins.
sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from decision_values import REPORTING_BINNING  # noqa: E402
from scale_cuts import (  # noqa: E402
    add_scale_cut_convention_argument,
    resolve_fiducial_bin_window,
    resolve_fiducial_scale_cuts,
)


BOOTSTRAP_SAMPLES = 10_000
BOOTSTRAP_SEED = 42


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-mock-statistics", type=Path, required=True)
    parser.add_argument("--pure-eb-pte", type=Path)
    parser.add_argument("--cosebis-pte", type=Path)
    parser.add_argument("--cl-pte", type=Path)
    parser.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)
    add_scale_cut_convention_argument(parser)
    parser.add_argument(
        "--observed-ptes", help="Optional JSON mapping statistic name to PTE override"
    )
    parser.add_argument("--global-alpha", type=float, default=0.05)
    parser.add_argument(
        "--pte-threshold",
        help="Adopted global PTE threshold (option ID or value), marked on the figure",
    )
    parser.add_argument(
        "--emit",
        choices=("calibration", "figure"),
        default="calibration",
        help=(
            "What this invocation writes: the calibration json (default), or "
            "the cumulative distribution of the per-mock minimum PTE."
        ),
    )
    parser.add_argument("--out-file", type=Path, required=True)
    return parser.parse_args()


def _wilson_68(count: int, n: int) -> list[float]:
    """Return the two-sided 68% Wilson interval for a binomial fraction."""
    if n <= 0:
        raise ValueError("cannot calculate an interval with no mocks")
    z = NormalDist().inv_cdf(0.84)
    fraction = count / n
    denominator = 1.0 + z * z / n
    centre = (fraction + z * z / (2.0 * n)) / denominator
    half_width = z * math.sqrt(
        fraction * (1.0 - fraction) / n + z * z / (4.0 * n * n)
    ) / denominator
    return [max(0.0, centre - half_width), min(1.0, centre + half_width)]


def _threshold_with_interval(values: np.ndarray, alpha: float) -> tuple[float, list[float]]:
    threshold = float(np.quantile(values, alpha))
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, values.size, size=(BOOTSTRAP_SAMPLES, values.size))
    bootstrap_thresholds = np.quantile(values[indices], alpha, axis=1)
    interval = [
        float(np.quantile(bootstrap_thresholds, 0.16)),
        float(np.quantile(bootstrap_thresholds, 0.84)),
    ]
    return threshold, interval


def _observed_ptes_from_artifacts(
    args: argparse.Namespace,
    stat_names: list[str],
    window: tuple[int, int],
    *,
    full_range: bool = False,
) -> dict[str, float]:
    missing = [name for name in ("pure_eb_pte", "cosebis_pte", "cl_pte") if getattr(args, name) is None]
    if missing:
        raise ValueError("artifact inputs required unless --observed-ptes is supplied: " + ", ".join(missing))
    window_min_bin, window_max_bin = window

    with np.load(args.pure_eb_pte, allow_pickle=False) as data:
        # These matrices are indexed [min_bin, max_bin] inclusively.  The
        # full-range row is the [0, 19] cell; the fiducial cell is resolved
        # from the declared scale-cut convention.
        pure_window = (0, -1) if full_range else (window_min_bin, window_max_bin)
        pure_values = {
            "xi+^B": float(data["pte_xip_B"][pure_window]),
            "xi-^B": float(data["pte_xim_B"][pure_window]),
            "xi_tot^B": float(data["pte_combined"][pure_window]),
        }

    with np.load(args.cosebis_pte, allow_pickle=False) as data:
        # The COSEBIS matrix uses theta-grid edges (21 entries), so the
        # reporting-bin window maps to the edge pair [min, max+1].
        cosebis_min = 0 if full_range else window_min_bin
        cosebis_max = -1 if full_range else window_max_bin + 1
        cosebis_values = {
            "B_n<=6": float(data["pte_B_6"][cosebis_min, cosebis_max]),
            "B_n<=20": float(data["pte_B_20"][cosebis_min, cosebis_max]),
        }

    with np.load(args.cl_pte, allow_pickle=False) as data:
        version = str(np.asarray(data["fiducial_version"]).item())
        ell = np.asarray(data[f"{version}__ell"], dtype=float)
        matrix = np.asarray(data[f"{version}__pte_matrix"], dtype=float)
        if full_range:
            lower, upper = 0, -1
        else:
            # For the adopted strict 300 < ell < 1600 cut, take [11, 27]:
            # ell[11]=327.5 is the first available lower edge above 300 and
            # ell[27]=1544 is the last available upper edge below 1600.
            lower = int(np.flatnonzero(ell > 300)[0])
            upper = int(np.flatnonzero(ell < 1600)[-1])
        cl_values = {"C_ell^BB": float(matrix[lower, upper])}

    observed = {**cosebis_values, **pure_values, **cl_values}
    if set(observed) != set(stat_names):
        raise ValueError("artifact-derived PTE names do not match stat_names")
    return {name: observed[name] for name in stat_names}


def _plot_min_pte(
    out_path: Path,
    one_sided: np.ndarray,
    observed_min: float,
    p_global: float,
    global_alpha: float,
    threshold: float,
    adopted_threshold: float,
    n_stats: int,
) -> None:
    """Cumulative distribution of the per-mock minimum PTE, read at 5%, at the adopted threshold and at the data."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    try:
        from plotting_utils import COLOR_B, FIG_WIDTH_SINGLE, PAPER_MPLSTYLE  # type: ignore

        plt.style.use(PAPER_MPLSTYLE)
    except Exception:  # pragma: no cover - style is cosmetic
        COLOR_B = "crimson"
        FIG_WIDTH_SINGLE = 3.54

    # The unbinned empirical CDF in percent: each rise is one realisation, and
    # its height at x is the suite's false-alarm rate for a threshold x.
    # Log-log, so the trend shows over several decades up to ~100%; x starts at
    # 1e-3, where the first few single-realisation steps have gone.
    n = one_sided.size
    x = np.concatenate([np.sort(one_sided), [1.0]])
    y = 100.0 * np.concatenate([np.arange(1, n + 1) / n, [1.0]])

    def cdf(value: float) -> float:
        return 100.0 * np.count_nonzero(one_sided <= value) / n

    x_lim = (1e-3, 0.3)
    y_lim = (0.85 * cdf(x_lim[0]), 100.0)

    # Drawn at the A&A column width, the width the paper prints it at, so the
    # style's font sizes are the ones the reader gets.
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_SINGLE, 2.8))
    ax.plot(x, y, drawstyle="steps-post", color="0.15", lw=1.1, zorder=2)
    # The expectation for n_stats independent tests, 1 - (1 - x)^n_stats.
    grid = np.geomspace(*x_lim, 400)
    ax.plot(
        grid, 100.0 * (1.0 - (1.0 - grid) ** n_stats),
        color="0.6", lw=0.8, ls=(0, (1, 1.5)), zorder=1.5,
    )
    # The curve read as drop lines: grey dashed from global_alpha across to the
    # calibrated threshold; solid grey from the adopted threshold up to the
    # curve and across to its false-alarm rate; crimson (the B mode's colour throughout the paper) from the data's
    # minimum up to p_global.
    ax.plot(
        [x_lim[0], threshold, threshold],
        [100.0 * global_alpha, 100.0 * global_alpha, y_lim[0]],
        color="0.55", lw=0.8, ls=(0, (4, 3)), zorder=1,
    )
    ax.plot(
        [adopted_threshold, adopted_threshold, x_lim[0]],
        [y_lim[0], cdf(adopted_threshold), cdf(adopted_threshold)],
        color="0.3", lw=0.7, zorder=1.8,
    )
    ax.plot(
        [observed_min, observed_min, x_lim[0]],
        [y_lim[0], 100.0 * p_global, 100.0 * p_global],
        color=COLOR_B, lw=1.0, zorder=3,
    )
    # Annotations at the axis-label size.
    size = plt.rcParams["axes.labelsize"]
    # "data" centred on the red horizontal leg (log-x midpoint), just above it.
    ax.text(
        np.sqrt(x_lim[0] * observed_min), 100.0 * p_global * 1.07, "data",
        color=COLOR_B, fontsize=size, ha="center", va="bottom",
    )
    ax.text(
        adopted_threshold * 1.08, y_lim[0] * 1.25, r"$p_\mathrm{thr}$",
        color="0.3", fontsize=size, ha="left", va="bottom",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.yaxis.set_major_locator(mticker.FixedLocator([1, 5, 10, 50, 100]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    # Open axes: left and bottom spines only, ticks only on those edges.
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(which="both", top=False, right=False)
    ax.set_xlabel(r"$\mathrm{PTE}_{\min}$")
    ax.set_ylabel(r"cumulative fraction of mocks (\%)")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}", flush=True)


def _calibration_payload(
    ptes: np.ndarray,
    stat_names: list[str],
    observed_ptes: dict[str, float],
    window: tuple[int, int],
    cut_arcmin: tuple[float, float],
    convention: str,
    global_alpha: float,
) -> dict:
    """Return one window's complete min-PTE calibration block."""

    one_sided = np.min(ptes, axis=1)
    two_sided = np.min(np.minimum(ptes, 1.0 - ptes), axis=1)
    observed_values = np.array([observed_ptes[name] for name in stat_names])
    observed_min = float(np.min(observed_values))
    observed_min_two_sided = float(
        np.min(np.minimum(observed_values, 1.0 - observed_values))
    )

    n_mocks = int(ptes.shape[0])
    one_count = int(np.count_nonzero(one_sided <= observed_min))
    two_count = int(np.count_nonzero(two_sided <= observed_min_two_sided))
    one_threshold, one_threshold_interval = _threshold_with_interval(
        one_sided, global_alpha
    )
    two_threshold, two_threshold_interval = _threshold_with_interval(
        two_sided, global_alpha
    )
    return {
        "stat_names": stat_names,
        "n_mocks": n_mocks,
        "observed_ptes": observed_ptes,
        "window_bins": list(window),
        "scale_cut_arcmin": list(cut_arcmin),
        "scale_cut_convention": convention,
        "observed_min": observed_min,
        "observed_min_two_sided": observed_min_two_sided,
        "p_global_one_sided": one_count / n_mocks,
        "p_global_one_sided_count": one_count,
        "p_global_one_sided_interval_68": _wilson_68(one_count, n_mocks),
        "p_global_two_sided": two_count / n_mocks,
        "p_global_two_sided_count": two_count,
        "p_global_two_sided_interval_68": _wilson_68(two_count, n_mocks),
        "global_alpha": global_alpha,
        "threshold_one_sided": one_threshold,
        "threshold_one_sided_interval_68": one_threshold_interval,
        "threshold_two_sided": two_threshold,
        "threshold_two_sided_interval_68": two_threshold_interval,
        "definition": (
            "Empirical fraction of clean mocks whose minimum PTE across the six "
            "correlated tests is <= the observed minimum."
        ),
    }


def main() -> None:
    args = _parse_args()
    if not 0.0 <= args.global_alpha <= 1.0:
        raise ValueError("--global-alpha must be between 0 and 1")

    with np.load(args.per_mock_statistics, allow_pickle=False) as data:
        ptes = np.asarray(data["pte"], dtype=float)
        ptes_full = np.asarray(data["pte_full"], dtype=float)
        stat_names = [str(name) for name in data["stat_names"]]

    for name, values in (("pte", ptes), ("pte_full", ptes_full)):
        if values.ndim != 2 or values.shape[1] != len(stat_names):
            raise ValueError(f"{name} must be a 2-D matrix with one column per stat_name")
        if values.shape[0] != ptes.shape[0]:
            raise ValueError("pte and pte_full must have the same number of mocks")
        if not np.all(np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
            raise ValueError(f"{name} values must be finite and lie in [0, 1]")

    window = resolve_fiducial_bin_window(args.reporting_binning, args.scale_cut_convention)
    cut_arcmin = resolve_fiducial_scale_cuts(args.reporting_binning, args.scale_cut_convention)
    full_window = (0, 19)
    full_cut_arcmin = args.reporting_binning[:2]

    if args.observed_ptes is None:
        observed_ptes = _observed_ptes_from_artifacts(args, stat_names, window)
        observed_ptes_full = _observed_ptes_from_artifacts(
            args, stat_names, window, full_range=True
        )
    else:
        try:
            observed_ptes = json.loads(args.observed_ptes)
        except json.JSONDecodeError as error:
            raise ValueError("--observed-ptes must be a JSON object") from error
        if not isinstance(observed_ptes, dict):
            raise ValueError("--observed-ptes must be a JSON object")
        if set(observed_ptes) != set(stat_names):
            raise ValueError("observed PTE names must match stat_names exactly")
        observed_ptes = {name: float(observed_ptes[name]) for name in stat_names}
        # There is no full-range override flag: artifact-derived values remain
        # the source for the second window even when the fiducial values are
        # supplied explicitly for a figure/debug invocation.
        observed_ptes_full = _observed_ptes_from_artifacts(
            args, stat_names, window, full_range=True
        )
    for label, values in (("observed_ptes", observed_ptes), ("observed_ptes_full", observed_ptes_full)):
        if not all(np.isfinite(value) and 0.0 <= value <= 1.0 for value in values.values()):
            raise ValueError(f"{label} values must be finite and lie in [0, 1]")

    fiducial = _calibration_payload(
        ptes,
        stat_names,
        observed_ptes,
        window,
        cut_arcmin,
        args.scale_cut_convention,
        args.global_alpha,
    )
    full_range = _calibration_payload(
        ptes_full,
        stat_names,
        observed_ptes_full,
        full_window,
        full_cut_arcmin,
        args.scale_cut_convention,
        args.global_alpha,
    )
    if args.emit == "figure":
        from numerics import resolve_pte_threshold  # type: ignore

        if args.pte_threshold is None:
            raise ValueError("--emit figure needs --pte-threshold")
        _plot_min_pte(
            args.out_file,
            np.min(ptes, axis=1),
            fiducial["observed_min"],
            fiducial["p_global_one_sided"],
            fiducial["global_alpha"],
            fiducial["threshold_one_sided"],
            resolve_pte_threshold(args.pte_threshold),
            len(stat_names),
        )
        return
    fiducial["full_range"] = full_range
    args.out_file.write_text(json.dumps(fiducial, indent=2) + "\n")


if __name__ == "__main__":
    main()
