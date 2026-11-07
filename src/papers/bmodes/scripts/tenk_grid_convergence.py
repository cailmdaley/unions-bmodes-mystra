"""Fine-grid convergence of the ADOPTED estimators: 1,000 vs 10,000 integration bins.

Sample variance cancels by construction: both grids are TreeCorr runs over the
same fiducial catalog on the same [0.5', 300'] extent, so the comparison is a
pure numerical-convergence statement about the estimators, not a statistical one.

What is compared
----------------
* **Pure modes** ξ_±^{E,B,amb} through the baseline chain — the fixed
  linear-transform operator evaluated on the fine grid, then collapsed to the
  20 reporting bins under the ``pure_eb_rebinning_weights`` rule.  The 10,000-bin
  vector is built by the very function the baseline output uses
  (``gather_pure_eb_chunks.data_vector``), so the only difference between the
  two sides is the grid.  The published pointwise ``dqags`` estimator is not
  exercised here.
* **COSEBIs** B_n, n ≤ ``nmodes``, at both scale cuts, with the mode integrals
  run on each grid.

Both grids are scored against the **same** covariance: the materialised 1,000-bin
baseline (the MC-propagated pure-mode covariance, and the T^T C_ξ T mode
covariance).  A 10,000-bin covariance does not exist and is not the question —
the question is whether the *data vector* moves.

Two single-file modes, mirroring the other paper scripts:
  compute:  --out-file <json>
  render:   --plot-from-json <json> --out-file <png>
"""

import argparse
import json

import numpy as np
from scipy import stats

_HAS_MPL = True
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from plotting_utils import FIG_WIDTH_SINGLE, PAPER_MPLSTYLE
except Exception:  # pragma: no cover - compute mode does not need a display stack
    _HAS_MPL = False

from gather_pure_eb_chunks import _load_xi, data_vector
from pure_mode_covariance import hartlap_samples
from pure_eb_transform import (
    add_integration_limits_argument,
    add_integration_method_argument,
    add_rebinning_weights_argument,
)
from scale_cuts import (
    add_scale_cut_convention_argument,
    reporting_grid_edges,
    resolve_fiducial_bin_window,
    resolve_fiducial_scale_cuts,
)
from catalog_versions import catalogs_from
from decision_values import FIDUCIAL_CATALOG, INTEGRATION_GRID, MODE_COUNT, REPORTING_BINNING

#: Pure-mode blocks, in the order the 6-block MC covariance stacks them.
EB_KEYS = ("xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")

#: The blocks this convergence test reports per bin.  E modes carry the signal
#: and are not a null test; the three below are what the paper makes claims about.
REPORTED_BLOCKS = ("xip_B", "xim_B", "xip_amb", "xim_amb")


def _chi2_pte(data, cov, n_eff):
    """χ² and its PTE, Hartlap-corrected when the covariance is MC-estimated.

    *n_eff* is the number of independent realisations behind the covariance, or
    ``None`` for a theory covariance, where no sample bias enters and the factor
    is 1 (``calculate_cosebis`` does the same).
    """
    dof = len(data)
    hartlap = 1.0 if n_eff is None else (n_eff - dof - 2) / (n_eff - 1)
    chi2 = float(hartlap * (data @ np.linalg.solve(cov, data)))
    return chi2, float(stats.chi2.sf(chi2, dof)), dof


def _fine_grid(path, nbins, min_sep, max_sep):
    """``_load_xi`` plus the log bin edges TreeCorr used for that grid."""
    gg = _load_xi(path, nbins)
    edges = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins + 1)
    gg["left_edges"] = edges[:-1]
    gg["right_edges"] = edges[1:]
    return gg


def _npairs_per_reporting_bin(gg_fine, reporting_edges):
    """Sum the fine grid's pair counts into the reporting bins.

    Uses the fine *node* separations, which is how the rebinning matrix assigns
    them.  A fine grid that resolves the reporting edges identically puts the
    same pairs in the same reporting bin whatever its resolution, so the two
    grids' sums agreeing is the check that both measurements are the same
    catalog over the same extent.
    """
    idx = np.digitize(gg_fine["meanr"], reporting_edges) - 1
    n_rep = len(reporting_edges) - 1
    sums = np.zeros(n_rep)
    counts = np.zeros(n_rep, dtype=int)
    inside = (idx >= 0) & (idx < n_rep)
    np.add.at(sums, idx[inside], gg_fine["npairs"][inside])
    np.add.at(counts, idx[inside], 1)
    return sums, counts


def _cosebis_Bn(gg_fine, nmodes, scale_cut, match_interval=None):
    """B_n (and E_n) from one fine grid at one scale cut.

    Same call chain as ``sp_validation.b_modes.calculate_cosebis``: conservative
    edge-based bin selection, then the COSEBI interval taken from the *retained
    bins' mean separations*.  That interval is therefore grid-dependent — the
    conservative selection lands on the first and last fine bin fully inside the
    cut, and where those fall moves with the resolution.

    COSEBIs are an orthogonal basis **on** [θ_min, θ_max], so two grids that
    retain different intervals are not two estimates of one number: they are two
    different statistics.  *match_interval* forces the selection to a given
    (θ_min, θ_max) by retaining the nodes inside it, which is what isolates the
    quadrature-resolution effect from the interval effect.
    """
    from cosmo_numba.B_modes.cosebis import COSEBIS

    if match_interval is not None:
        lo, hi = match_interval
        inds = np.flatnonzero(
            (gg_fine["meanr"] >= lo * (1 - 1e-12))
            & (gg_fine["meanr"] <= hi * (1 + 1e-12))
        )
        start, stop = int(inds[0]), int(inds[-1] + 1)
        theta = gg_fine["meanr"][inds]
        cosebis = COSEBIS(
            theta_min=float(theta.min()),
            theta_max=float(theta.max()),
            N_max=nmodes,
            precision=120,
        )
        En, Bn = cosebis.cosebis_from_xipm(
            theta, gg_fine["xip"][inds], gg_fine["xim"][inds], parallel=True
        )
        return (
            np.asarray(En),
            np.asarray(Bn),
            (float(theta.min()), float(theta.max())),
            (start, stop),
        )

    lo, hi = scale_cut
    start = int(np.searchsorted(gg_fine["left_edges"], lo, side="left"))
    stop = int(np.searchsorted(gg_fine["right_edges"], hi, side="right"))
    inds = np.arange(start, stop)
    theta = gg_fine["meanr"][inds]
    cosebis = COSEBIS(
        theta_min=float(np.min(theta)),
        theta_max=float(np.max(theta)),
        N_max=nmodes,
        precision=120,
    )
    En, Bn = cosebis.cosebis_from_xipm(
        theta, gg_fine["xip"][inds], gg_fine["xim"][inds], parallel=True
    )
    return np.asarray(En), np.asarray(Bn), (float(theta.min()), float(theta.max())), (start, stop)


def compute(
    version,
    nmodes,
    integration,
    reporting,
    *,
    xi_reporting,
    xi_1k,
    xi_10k,
    pure_eb_baseline,
    cosebis_baseline,
    out_file,
    scale_cut_convention,
    integration_method,
    rebinning_weights,
    integration_limits,
):
    min_sep, max_sep, nbins = reporting
    min_sep_int, max_sep_int, nbins_1k = integration
    nbins_10k = INTEGRATION_GRID["bins_10000"].nbins

    gg = _load_xi(xi_reporting, nbins)
    gg_1k = _fine_grid(xi_1k, nbins_1k, min_sep_int, max_sep_int)
    gg_10k = _fine_grid(xi_10k, nbins_10k, min_sep_int, max_sep_int)

    payload = {
        "meta": {
            "version": version,
            "nbins_reporting": nbins,
            "nbins_fine": {"coarse": nbins_1k, "fine": nbins_10k},
            "fine_extent_arcmin": [min_sep_int, max_sep_int],
            "decisions": {
                "scale_cut_convention": scale_cut_convention,
                "pure_eb_integration_method": integration_method,
                "pure_eb_rebinning_weights": rebinning_weights,
                "pure_eb_integration_limits": integration_limits,
            },
            "covariance": (
                "the 1,000-bin baseline is used for BOTH grids: the MC-propagated "
                "6-block pure-mode covariance, and the COSEBI T^T C_xi T mode "
                "covariance. No 10,000-bin covariance exists."
            ),
        }
    }

    # ---- grid sanity -------------------------------------------------------
    reporting_edges = reporting_grid_edges(reporting)
    npairs_1k, nodes_1k = _npairs_per_reporting_bin(gg_1k, reporting_edges)
    npairs_10k, nodes_10k = _npairs_per_reporting_bin(gg_10k, reporting_edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        npairs_rel = np.where(npairs_1k > 0, npairs_10k / npairs_1k - 1.0, np.nan)
    # How sharply does each grid resolve a reporting edge? Neither fine grid's
    # edges coincide with the reporting edges, so exactly one fine bin straddles
    # each of them, and that bin is assigned wholly to whichever reporting bin
    # its node falls in. The misassignment is that bin's pair count as a
    # fraction of the reporting bin's, and it shrinks with resolution — it is
    # the leading reason the two grids' per-bin pair sums differ at all. The
    # window itself never depends on this: it is resolved from the reporting
    # EDGES, never from bin means.
    edge_resolution = {}
    for label, gg_f in (("coarse", gg_1k), ("fine", gg_10k)):
        left, right = gg_f["left_edges"], gg_f["right_edges"]
        counts, fracs = [], []
        for e in reporting_edges:
            hit = np.flatnonzero((left < e) & (right > e))
            counts.append(int(hit.size))
            if hit.size:
                node = gg_f["meanr"][hit[0]]
                rep = int(np.digitize(node, reporting_edges) - 1)
                total = (
                    npairs_1k[rep] if label == "coarse" else npairs_10k[rep]
                ) if 0 <= rep < len(reporting_edges) - 1 else np.nan
                fracs.append(float(gg_f["npairs"][hit[0]] / total))
            else:
                fracs.append(0.0)
        edge_resolution[label] = {
            "fine_bins_straddling_each_reporting_edge": counts,
            "straddling_bin_npairs_fraction_of_its_reporting_bin": fracs,
        }
    first_bin, last_bin = resolve_fiducial_bin_window(reporting, scale_cut_convention)
    payload["grid_checks"] = {
        "reporting_edges_arcmin": reporting_edges.tolist(),
        "fiducial_window_bins_inclusive": [int(first_bin), int(last_bin)],
        "fiducial_window_theta_arcmin": [
            float(reporting_edges[first_bin]),
            float(reporting_edges[last_bin + 1]),
        ],
        "npairs_per_reporting_bin": {
            "coarse": npairs_1k.tolist(),
            "fine": npairs_10k.tolist(),
            "fine_over_coarse_minus_one": npairs_rel.tolist(),
        },
        "fine_nodes_per_reporting_bin": {
            "coarse": nodes_1k.tolist(),
            "fine": nodes_10k.tolist(),
        },
        "npairs_total": {
            "coarse": float(gg_1k["npairs"].sum()),
            "fine": float(gg_10k["npairs"].sum()),
            "fine_over_coarse_minus_one": float(
                gg_10k["npairs"].sum() / gg_1k["npairs"].sum() - 1.0
            ),
        },
        "reporting_edge_resolution": edge_resolution,
    }

    # ---- pure modes --------------------------------------------------------
    base = np.load(pure_eb_baseline)
    cov_pure = base["cov_pure_eb"]
    dv_kwargs = dict(
        min_sep=min_sep,
        max_sep=max_sep,
        nbins=nbins,
        integration_method=integration_method,
        rebinning_weights=rebinning_weights,
        integration_limits=integration_limits,
    )
    print("Pure modes on the 1,000-bin grid (recompute, checked against baseline)...")
    modes_1k = data_vector(gg, gg_1k, **dv_kwargs)
    print("Pure modes on the 10,000-bin grid...")
    modes_10k = data_vector(gg, gg_10k, **dv_kwargs)

    recompute_check = {}
    for key, values in zip(EB_KEYS, modes_1k):
        stored = base[key]
        recompute_check[key] = float(np.max(np.abs(values - stored)))
    payload["meta"]["coarse_recompute_max_abs_diff_vs_baseline"] = recompute_check
    scale = float(np.max(np.abs(base["xip_total"])))
    worst = max(recompute_check.values())
    if worst > 1e-10 * scale:
        raise SystemExit(
            "the 1,000-bin recompute does not reproduce the materialised baseline "
            f"pure-mode vector (max |diff| = {worst:.3e}, max|xi_+| = {scale:.3e}); "
            "the comparison would not be grid-only."
        )

    pure = {}
    sigma_by_key = {}
    for i, key in enumerate(EB_KEYS):
        block = cov_pure[i * nbins : (i + 1) * nbins, i * nbins : (i + 1) * nbins]
        sigma_by_key[key] = np.sqrt(np.diag(block))
    for key in REPORTED_BLOCKS:
        i = EB_KEYS.index(key)
        v1, v10 = np.asarray(modes_1k[i]), np.asarray(modes_10k[i])
        sigma = sigma_by_key[key]
        pure[key] = {
            "theta_arcmin": gg["meanr"].tolist(),
            "sigma": sigma.tolist(),
            "value_coarse": v1.tolist(),
            "value_fine": v10.tolist(),
            "diff": (v10 - v1).tolist(),
            "diff_over_sigma": ((v10 - v1) / sigma).tolist(),
        }
    payload["pure_modes"] = pure

    # Fiducial-window PTEs, both grids, the 1,000-bin covariance for both, with
    # the Hartlap factor that covariance artifact records (None: exact).
    n_eff_mc = hartlap_samples(base)
    win = slice(first_bin, last_bin + 1)
    pure_ptes = {}
    for key in ("xip_B", "xim_B"):
        i = EB_KEYS.index(key)
        block = cov_pure[i * nbins : (i + 1) * nbins, i * nbins : (i + 1) * nbins][
            win, win
        ]
        entry = {}
        for label, modes in (("coarse", modes_1k), ("fine", modes_10k)):
            chi2, pte, dof = _chi2_pte(np.asarray(modes[i])[win], block, n_eff_mc)
            entry[label] = {"chi2": chi2, "pte": pte, "dof": dof}
        pure_ptes[key] = entry
    # Joint ξ_+^B ⊕ ξ_-^B, with the cross block, as calculate_eb_statistics does.
    ip, im = EB_KEYS.index("xip_B"), EB_KEYS.index("xim_B")
    sp = slice(ip * nbins + first_bin, ip * nbins + last_bin + 1)
    sm = slice(im * nbins + first_bin, im * nbins + last_bin + 1)
    cov_joint = np.block(
        [
            [cov_pure[sp, sp], cov_pure[sp, sm]],
            [cov_pure[sm, sp], cov_pure[sm, sm]],
        ]
    )
    entry = {}
    for label, modes in (("coarse", modes_1k), ("fine", modes_10k)):
        joint = np.concatenate(
            [np.asarray(modes[ip])[win], np.asarray(modes[im])[win]]
        )
        chi2, pte, dof = _chi2_pte(joint, cov_joint, n_eff_mc)
        entry[label] = {"chi2": chi2, "pte": pte, "dof": dof}
    pure_ptes["combined"] = entry

    # ---- COSEBIs -----------------------------------------------------------
    cos_base = np.load(cosebis_baseline, allow_pickle=True)
    scale_cuts = {
        "fiducial": resolve_fiducial_scale_cuts(reporting, scale_cut_convention),
        "full": (
            reporting.min_sep,
            reporting.max_sep,
        ),
    }
    cosebis_payload = {}
    cosebis_ptes = {}
    for key, cut in scale_cuts.items():
        stored_cut = tuple(float(v) for v in cos_base[f"{key}_scale_cut"])
        if not np.allclose(stored_cut, cut):
            raise SystemExit(
                f"{key} scale cut {cut} does not match the baseline COSEBI NPZ "
                f"{stored_cut}; the two sides would not be comparable."
            )
        Bn_1k = np.asarray(cos_base[f"{key}_Bn"])
        cov_B = np.asarray(cos_base[f"{key}_cov"])[nmodes:, nmodes:]
        sigma_B = np.sqrt(np.diag(cov_B))
        print(f"COSEBIs on the 10,000-bin grid, {key} cut {cut}...")
        _, Bn_10k, interval_10k, bins_10k = _cosebis_Bn(gg_10k, nmodes, cut)
        # The 1,000-bin interval, for the record; B_n itself is the baseline's.
        start = int(np.searchsorted(gg_1k["left_edges"], cut[0], side="left"))
        stop = int(np.searchsorted(gg_1k["right_edges"], cut[1], side="right"))
        theta_1k = gg_1k["meanr"][start:stop]
        interval_1k = (float(theta_1k.min()), float(theta_1k.max()))
        # The same 10,000-bin ξ± integrated over the 1,000-bin grid's OWN
        # retained interval. COSEBIs are an orthogonal basis on [θ_min, θ_max],
        # so `fine` above differs from `coarse` by resolution AND by interval;
        # this third vector holds the interval fixed, leaving resolution alone.
        print(f"COSEBIs on the 10,000-bin grid, {key}, interval matched to 1k...")
        _, Bn_10k_matched, interval_matched, bins_matched = _cosebis_Bn(
            gg_10k, nmodes, cut, match_interval=interval_1k
        )
        # The fourth corner: the 1,000-bin ξ± over the 10,000-bin grid's
        # interval. The coarse grid cannot land on [12.003', 82.800'] exactly —
        # its nodes are 0.64% apart — so this retains the coarse nodes that fall
        # inside that interval, which reaches further up than the coarse grid's
        # own conservative selection does. The achieved interval is recorded.
        print(f"COSEBIs on the 1,000-bin grid, {key}, interval matched to 10k...")
        _, Bn_1k_matched, interval_1k_matched, bins_1k_matched = _cosebis_Bn(
            gg_1k, nmodes, cut, match_interval=interval_10k
        )
        Bn_10k = np.asarray(Bn_10k)
        Bn_10k_matched = np.asarray(Bn_10k_matched)
        Bn_1k_matched = np.asarray(Bn_1k_matched)
        cosebis_payload[key] = {
            "scale_cut_arcmin": list(cut),
            "n": list(range(1, nmodes + 1)),
            "sigma": sigma_B.tolist(),
            "value_coarse": Bn_1k.tolist(),
            "value_fine": Bn_10k.tolist(),
            "value_fine_interval_matched": Bn_10k_matched.tolist(),
            "value_coarse_interval_matched": Bn_1k_matched.tolist(),
            "diff_over_sigma_coarse_interval_matched": (
                (Bn_1k_matched - Bn_1k) / sigma_B
            ).tolist(),
            "diff": (Bn_10k - Bn_1k).tolist(),
            "diff_over_sigma": ((Bn_10k - Bn_1k) / sigma_B).tolist(),
            "diff_over_sigma_interval_matched": (
                (Bn_10k_matched - Bn_1k) / sigma_B
            ).tolist(),
            "diff_over_sigma_interval_only": (
                (Bn_10k - Bn_10k_matched) / sigma_B
            ).tolist(),
            "effective_interval_arcmin": {
                "coarse": list(interval_1k),
                "fine": list(interval_10k),
                "fine_interval_matched": list(interval_matched),
                "coarse_interval_matched": list(interval_1k_matched),
            },
            "retained_fine_bins": {
                "coarse": [start, stop],
                "fine": [int(bins_10k[0]), int(bins_10k[1])],
                "fine_interval_matched": [
                    int(bins_matched[0]),
                    int(bins_matched[1]),
                ],
                "coarse_interval_matched": [
                    int(bins_1k_matched[0]),
                    int(bins_1k_matched[1]),
                ],
            },
        }
        entry = {}
        for n_sub in (6, nmodes):
            sub = {}
            for label, Bn in (
                ("coarse", Bn_1k),
                ("coarse_interval_matched", Bn_1k_matched),
                ("fine", Bn_10k),
                ("fine_interval_matched", Bn_10k_matched),
            ):
                chi2, pte, dof = _chi2_pte(
                    Bn[:n_sub], cov_B[:n_sub, :n_sub], None
                )
                sub[label] = {"chi2": chi2, "pte": pte, "dof": dof}
            entry[f"n_le_{n_sub}"] = sub
        cosebis_ptes[key] = entry
    payload["cosebis"] = cosebis_payload

    payload["ptes"] = {
        "pure_modes_fiducial_window": pure_ptes,
        "cosebis": cosebis_ptes,
    }

    with open(out_file, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Saved {out_file}")


#: The four grid x interval corners, in plot order, as (payload key, legend).
COSEBI_CURVES = (
    ("coarse", "1,000 bins, own interval"),
    ("fine", "10,000 bins, own interval"),
    ("fine_interval_matched", "10,000 bins, 1,000-bin interval"),
    ("coarse_interval_matched", "1,000 bins, 10,000-bin interval"),
)


def render_cosebis(json_path, out_file):
    """B_n / σ_n itself, for all four grid × interval corners, both scale cuts.

    The difference figure answers "how far does it move"; this one answers
    "what are the modes", which is what a reader comparing against the paper's
    data-vector figure needs. σ_n is the 1,000-bin T^T C_ξ T mode covariance
    throughout — the only one that exists — so the ±1 band is one and the same
    for every curve.
    """

    if not _HAS_MPL:  # pragma: no cover
        raise SystemExit("matplotlib/plotting_utils unavailable in this environment")
    with open(json_path) as f:
        payload = json.load(f)
    plt.style.use(PAPER_MPLSTYLE)

    panels = [
        ("fiducial", "fiducial"),
        ("full", "full"),
    ]
    colors = sns.color_palette("colorblind", 4)
    markers = ("o", "s", "^", "D")
    offsets = (-0.24, -0.08, 0.08, 0.24)

    fig, axes = plt.subplots(
        1, 2, figsize=(FIG_WIDTH_SINGLE * 2.3, FIG_WIDTH_SINGLE * 0.78), sharey=True
    )
    for ax, (key, _) in zip(axes, panels):
        cos = payload["cosebis"][key]
        n = np.asarray(cos["n"], dtype=float)
        sigma = np.asarray(cos["sigma"])
        intervals = cos["effective_interval_arcmin"]
        ax.axhspan(-1, 1, color="0.92", zorder=0)
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.6)
        ax.axvspan(0.5, 6.5, color="0.96", zorder=0)
        for (curve, legend), color, marker, dx in zip(
            COSEBI_CURVES, colors, markers, offsets
        ):
            values = np.asarray(cos[f"value_{curve}"]) / sigma
            ilo, ihi = intervals[curve]
            ax.errorbar(
                n + dx,
                values,
                yerr=np.ones_like(values),
                fmt=marker,
                color=color,
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=0.4,
                markersize=3.6,
                capsize=1.5,
                capthick=0.7,
                elinewidth=0.7,
                linewidth=0,
                label=rf"{legend} [${ilo:.3f}$, ${ihi:.3f}'$]",
            )
        # Headroom for the four-line legend, which would otherwise sit on the
        # low-n modes where the interesting structure is.
        stacked = np.concatenate(
            [np.asarray(cos[f"value_{c}"]) / sigma for c, _ in COSEBI_CURVES]
        )
        span = stacked.max() + 1 - (stacked.min() - 1)
        ax.set_ylim(stacked.min() - 1 - 0.08 * span, stacked.max() + 1 + 0.55 * span)
        cut = cos["scale_cut_arcmin"]
        ax.set_xlabel("COSEBI mode $n$")
        ax.set_xlim(0.4, len(n) + 0.6)
        ax.set_xticks(np.arange(1, len(n) + 1))
        ax.set_xticklabels(
            [str(i) for i in range(1, len(n) + 1)],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=6,
        )
        ax.set_title(rf"$\theta = {cut[0]:.0f}$--${cut[1]:.0f}'$", fontsize=9)
        ax.legend(loc="upper left", fontsize=5.4, frameon=True, framealpha=0.9)

    axes[0].set_ylabel(r"$B_n / \sigma_n$")
    fig.suptitle(
        r"COSEBI $B_n$ across integration grid $\times$ COSEBI interval. "
        r"$\sigma_n$ is the 1,000-bin $T^{\mathsf{T}} C_\xi T$ mode covariance "
        "for every curve.",
        fontsize=7.5,
    )
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_file}")


def render(json_path, out_file):
    if not _HAS_MPL:  # pragma: no cover
        raise SystemExit("matplotlib/plotting_utils unavailable in this environment")
    with open(json_path) as f:
        payload = json.load(f)
    plt.style.use(PAPER_MPLSTYLE)

    lo, hi = payload["grid_checks"]["fiducial_window_theta_arcmin"]
    labels = {
        "xip_B": r"$\xi_+^{B}$",
        "xim_B": r"$\xi_-^{B}$",
        "xip_amb": r"$\xi_+^{\rm amb}$",
        "xim_amb": r"$\xi_-^{\rm amb}$",
    }
    colors = sns.color_palette("colorblind", 5)

    fig, axes = plt.subplots(
        1, 5, figsize=(FIG_WIDTH_SINGLE * 2.6, FIG_WIDTH_SINGLE * 0.62)
    )
    for ax in axes[1:4]:
        ax.sharey(axes[0])
    for ax, (key, label), color in zip(axes[:4], labels.items(), colors):
        block = payload["pure_modes"][key]
        theta = np.asarray(block["theta_arcmin"])
        d = np.asarray(block["diff_over_sigma"])
        ax.axhspan(-0.1, 0.1, color="0.92", zorder=0)
        ax.axvspan(lo, hi, color=colors[4], alpha=0.13, zorder=0)
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.6)
        ax.plot(theta, d, "o-", color=color, markersize=3.5, linewidth=0.9)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\theta$ [arcmin]")
        ax.set_title(label, fontsize=9)

    ax = axes[4]
    cos = payload["cosebis"]["fiducial"]
    n = np.asarray(cos["n"])
    ax.axhspan(-0.1, 0.1, color="0.92", zorder=0)
    ax.axhline(0, color="black", linewidth=0.8, alpha=0.6)
    ax.plot(
        n,
        np.asarray(cos["diff_over_sigma"]),
        "s-",
        color=colors[3],
        markersize=3.5,
        linewidth=0.9,
        label="as measured",
    )
    ax.plot(
        n,
        np.asarray(cos["diff_over_sigma_interval_matched"]),
        "^--",
        color=colors[2],
        markersize=3.5,
        linewidth=0.9,
        label="interval matched",
    )
    ax.set_xlabel("COSEBI mode $n$")
    ax.set_title(rf"$B_n$ (${lo:.0f}$--${hi:.0f}'$)", fontsize=9)
    span = np.concatenate(
        [
            np.asarray(cos["diff_over_sigma"]),
            np.asarray(cos["diff_over_sigma_interval_matched"]),
        ]
    )
    pad = 0.45 * (span.max() - span.min())
    ax.set_ylim(span.min() - 0.1 * pad, span.max() + pad)
    ax.legend(loc="upper left", fontsize=6, frameon=True, framealpha=0.9)
    for a in axes[1:4]:
        a.tick_params(labelleft=False)

    axes[0].set_ylabel(r"$(10{,}000\ \mathrm{bins} - 1{,}000\ \mathrm{bins})\,/\,\sigma$")
    fig.suptitle(
        "Integration-grid convergence of the adopted estimators "
        rf"(shaded band: fiducial ${lo:.0f}$--${hi:.0f}'$ window; "
        r"grey: $\pm 0.1\,\sigma$)",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_file}")


def _from_cli(argv=None):

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fiducial-catalog", choices=FIDUCIAL_CATALOG)
    ap.add_argument("--leakage-correction", choices=["applied", "uncorrected"])
    ap.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__)
    ap.add_argument("--integration-grid", type=INTEGRATION_GRID.__getitem__)
    ap.add_argument("--mode-count", type=MODE_COUNT.__getitem__)
    ap.add_argument("--xi-reporting")
    ap.add_argument("--xi-1k", help="1,000-bin integration-grid TreeCorr dump")
    ap.add_argument("--xi-10k", help="10,000-bin integration-grid TreeCorr dump")
    ap.add_argument("--pure-eb-baseline", help="Baseline pure_eb_data NPZ")
    ap.add_argument("--cosebis-baseline", help="Baseline cosebis_modes_data NPZ")
    ap.add_argument("--plot-from-json", help="Render mode: the computed JSON")
    ap.add_argument(
        "--plot-cosebis-from-json",
        help="Render mode: B_n itself for the four grid x interval corners",
    )
    ap.add_argument("--out-file", required=True, help="lc {output}")
    add_scale_cut_convention_argument(ap)
    add_integration_method_argument(ap)
    add_integration_limits_argument(ap)
    add_rebinning_weights_argument(ap)
    a = ap.parse_args(argv)

    if a.plot_cosebis_from_json:
        render_cosebis(a.plot_cosebis_from_json, a.out_file)
        return
    if a.plot_from_json:
        render(a.plot_from_json, a.out_file)
        return
    required = ("fiducial_catalog", "leakage_correction", "reporting_binning", "integration_grid",
                "mode_count", "xi_reporting", "xi_1k", "xi_10k", "pure_eb_baseline",
                "cosebis_baseline")
    missing = [f"--{r.replace('_', '-')}" for r in required if not getattr(a, r)]
    if missing:
        ap.error(f"compute mode needs {', '.join(missing)}")
    compute(
        catalogs_from(a).fiducial,
        a.mode_count,
        a.integration_grid,
        a.reporting_binning,
        xi_reporting=a.xi_reporting,
        xi_1k=a.xi_1k,
        xi_10k=a.xi_10k,
        pure_eb_baseline=a.pure_eb_baseline,
        cosebis_baseline=a.cosebis_baseline,
        out_file=a.out_file,
        scale_cut_convention=a.scale_cut_convention,
        integration_method=a.integration_method,
        rebinning_weights=a.rebinning_weights,
        integration_limits=a.integration_limits,
    )


if __name__ == "__main__":
    _from_cli()
