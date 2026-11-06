"""The decision scan's one computed table: every variant's verdict, in one json.

Reads the scan's two upstream artifacts -- the per-variant covariance bundle
and the per-variant mock pure-mode ensemble -- and answers, per variant, the
three questions that decide between the variants:

* **mock bias.**  The mocks carry no B power by construction, so the
  across-realisation mean of xi_+^B / xi_-^B is that variant's own
  rectification.  Reported as chi2/PTE of the mean against the empirical
  covariance of the mean (sigma_emp/sqrt(N)), Hartlap-corrected, plus the
  per-bin mean in units of the analytic sigma.
* **calibration.**  sigma_emp/sigma_ana per bin, and the KS p-value of the
  across-realisation PTE distribution against uniform for xi_+^B, xi_-^B and
  xi_tot^B.  A variant can be unbiased and still mis-calibrated if its
  propagated covariance has the wrong amplitude.
* **the data.**  The real catalogue's PTE for the same three statistics at the
  adopted window, computed through that variant's own covariance -- so the
  mock verdict and the number the paper would report are read off one table.

Every per-mock chi2/PTE reuses ``mock_campaign_analysis``'s own covariance
slicing, scale-cut constants and Hartlap convention rather than restating
them, and the data PTE reuses ``calculate_pure_eb_ptes``; this stage adds the
variant loop and nothing else.

The Hartlap factor is NOT the pipeline's HARTLAP_N_PURE = 2000.  The scan's
covariances are built from 500 draws, so the correction on their inverse must
use 500; the count is read from each bundle's own provenance rather than
restated here, so the two can never drift apart.

    python pure_eb_scan_statistics.py --covariance-bundle B1.tar.gz B2.tar.gz ... \
        --mock-modes M.npz --cov-integration C --reference-covariance mockmatched \
        --out-file S.json
"""

from __future__ import annotations

import argparse
import json
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
from scipy import stats

from mock_campaign_analysis import (
    PURE_THETA_BY_CONVENTION,
    _apply_scale_cut_convention,
    _load_pure_covariance,
)
from plotting_utils import compute_chi2_pte

# The two script trees are separate sys.path roots; the data-side PTE
# calculation lives with the paper-facing scripts and is shared rather than
# duplicated, exactly as mock_pure_eb_scatter shares pure_eb_transform.
sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from calculate_pure_eb_ptes import calculate_ptes  # noqa: E402
from gather_pure_eb_chunks import _load_xi, data_vector  # noqa: E402
from decision_values import REPORTING_BINNING  # noqa: E402
from scale_cuts import resolve_fiducial_bin_window  # noqa: E402
from pure_eb_scan_variants import (  # noqa: E402
    ADOPTED,
    BY_NAME,
    DROPPED,
    PUBLISHED,
    add_variants_argument,
    resolve,
)

STATS = ("xi+^B", "xi-^B", "xi_tot^B")
# calculate_pure_eb_ptes' matrix names, indexed [min_bin, max_bin] inclusive.
DATA_PTE_KEYS = {"xi+^B": "pte_xip_B", "xi-^B": "pte_xim_B", "xi_tot^B": "pte_combined"}


#: The six pure-mode blocks a bundle's data npz carries, in reference order.
_EB_KEYS = ("xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")


def _rebuild_data_npz(variant, donor_npz, gg, gg_int, out_path, *, min_sep, max_sep, nbins):
    """Write *variant*'s own data vector against the donor's covariance.

    A variant with ``covariance_from`` set shares its donor's propagated
    covariance but not its data vector -- the data vector is exactly what such
    a variant changes.  So the borrowed ``<donor>__data.npz`` is copied with
    its six mode blocks replaced by this variant's own, computed through
    ``gather_pure_eb_chunks.data_vector``: the same function the bundle's own
    data vector came from, so the two differ only in the variant's parameters.

    Returns the recomputed six-tuple, so the caller can check the shared code
    path against the stored vector for a variant that has both.
    """

    modes = data_vector(
        gg,
        gg_int,
        min_sep=min_sep,
        max_sep=max_sep,
        nbins=nbins,
        integration_method=variant.integration_method,
        rebinning_weights=variant.rebinning_weights,
        integration_limits=variant.integration_limits,
        geometry=variant.geometry,
        local_from_eval_grid=variant.local_from_eval_grid,
    )
    if out_path is not None:
        with np.load(donor_npz, allow_pickle=False) as donor:
            package = {key: np.asarray(donor[key]) for key in donor.files}
        package.update(dict(zip(_EB_KEYS, (np.asarray(m, dtype=float) for m in modes))))
        np.savez(out_path, **package)
    return modes


def _finite(value):
    """JSON-safe float: NaN and +/-inf become None rather than invalid JSON."""
    value = float(value)
    return value if np.isfinite(value) else None


def _mean_vector_chi2(vectors: np.ndarray) -> dict[str, object]:
    """chi2/PTE of the mean against the empirical covariance of the mean."""

    n_real, n_bins = vectors.shape
    mean = vectors.mean(axis=0)
    cov_mean = np.cov(vectors, rowvar=False, ddof=1) / n_real
    # Hartlap on an N-realisation estimate of an n_bins x n_bins covariance.
    hartlap = (n_real - n_bins - 2) / (n_real - 1)
    chi2 = float(mean @ np.linalg.inv(cov_mean) @ mean * hartlap)
    return {
        "chi2": _finite(chi2),
        "dof": int(n_bins),
        "pte": _finite(stats.chi2.sf(chi2, n_bins)),
        "hartlap": _finite(hartlap),
        "n_realizations": int(n_real),
    }


def _variant_statistics(modes, variant, cov_path, hartlap_n):
    theta = np.asarray(modes[f"{variant.name}__theta"])[0]
    xip = np.asarray(modes[f"{variant.name}__xip_B"])
    xim = np.asarray(modes[f"{variant.name}__xim_B"])
    n_real = xip.shape[0]

    covariance = _load_pure_covariance(cov_path, theta)
    mask = covariance["mask"]
    cut = {
        "xi+^B": (xip[:, mask], covariance["xip_cut"]),
        "xi-^B": (xim[:, mask], covariance["xim_cut"]),
        "xi_tot^B": (
            np.concatenate([xip[:, mask], xim[:, mask]], axis=1),
            covariance["joint_cut"],
        ),
    }

    per_mock = {}
    for name, (vectors, block) in cut.items():
        values = [
            compute_chi2_pte(vector, block, n_samples=hartlap_n)
            for vector in vectors
        ]
        chi2 = np.array([entry[0] for entry in values], dtype=float)
        pte = np.array([entry[1] for entry in values], dtype=float)
        dof = int(values[0][2])
        finite = np.isfinite(chi2) & np.isfinite(pte)
        ks = stats.kstest(pte[finite], "uniform") if finite.any() else None
        per_mock[name] = {
            "dof": dof,
            "n": int(finite.sum()),
            "hartlap_n": int(hartlap_n),
            "chi2": [_finite(value) for value in chi2],
            "pte": [_finite(value) for value in pte],
            "chi2_mean": _finite(chi2[finite].mean()) if finite.any() else None,
            "expected_mean_nu": dof,
            "ks_pvalue": _finite(ks.pvalue) if ks else None,
            "ks_statistic": _finite(ks.statistic) if ks else None,
            "mean_vector": _mean_vector_chi2(vectors),
        }

    # Bias and calibration on the FULL reporting range, so the figure can shade
    # the cut rather than being restricted to it.
    profiles = {}
    for name, values, sigma_ana in (
        ("xi+^B", xip, covariance["sigma_xip_full"]),
        ("xi-^B", xim, covariance["sigma_xim_full"]),
    ):
        sigma_emp = values.std(axis=0, ddof=1)
        profiles[name] = {
            "theta": [_finite(value) for value in theta],
            "mean_over_sigma_ana": [
                _finite(value) for value in values.mean(axis=0) / sigma_ana
            ],
            # The error on that mean, in the same units.
            "stderr_over_sigma_ana": [
                _finite(value) for value in sigma_emp / sigma_ana / np.sqrt(n_real)
            ],
            "sigma_ratio": [_finite(value) for value in sigma_emp / sigma_ana],
            # 1/sqrt(2(N-1)) is the sampling error of an N-realisation stddev.
            "sigma_ratio_error": _finite(1.0 / np.sqrt(2 * (n_real - 1))),
        }

    return {"per_mock": per_mock, "profiles": profiles, "n_realizations": int(n_real)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--covariance-bundle",
        nargs="+",
        required=True,
        help=(
            "One per-variant covariance bundle from pure_eb_scan_covariances. "
            "Repeatable; every variant in --variants must be covered."
        ),
    )
    parser.add_argument("--mock-modes", required=True)
    parser.add_argument("--cov-integration", required=True)
    parser.add_argument(
        "--reference-covariance",
        choices=("fiducial", "mockmatched"),
        required=True,
        help=(
            "Which of the variant's two covariances the mock chi2/PTE run "
            "against: 'fiducial' is the data's own, 'mockmatched' the same "
            "propagation at the mocks' noise."
        ),
    )
    parser.add_argument(
        "--scale-cut-convention",
        choices=tuple(PURE_THETA_BY_CONVENTION),
        default="grid_bin_edges",
    )
    parser.add_argument("--reporting-binning", type=REPORTING_BINNING.__getitem__, required=True)
    parser.add_argument("--npatch", type=int, default=1)
    parser.add_argument(
        "--xi-reporting",
        required=True,
        help=(
            "Measured 20-bin reporting xi. Needed to rebuild the data vector "
            "of a variant whose covariance is aliased to another's, whose "
            "bundle therefore carries the donor's data vector."
        ),
    )
    parser.add_argument(
        "--xi-integration",
        required=True,
        help="Measured 1000-bin integration-grid xi; see --xi-reporting.",
    )
    parser.add_argument("--min-sep", type=float, default=1.0)
    parser.add_argument("--max-sep", type=float, default=250.0)
    parser.add_argument("--nbins", type=int, default=20)
    parser.add_argument("--nbins-int", type=int, default=1000)
    parser.add_argument("--out-file", required=True)
    add_variants_argument(parser)
    args = parser.parse_args()

    # One definition of the fiducial window, shared with the paper's PTE
    # matrices; see scale_cuts.resolve_fiducial_bin_window.
    window_min_bin, window_max_bin = resolve_fiducial_bin_window(
        args.reporting_binning, args.scale_cut_convention
    )
    # _load_pure_covariance masks with mock_campaign_analysis's module-level
    # theta cut, so the convention has to be bound before the first call.
    _apply_scale_cut_convention(args.scale_cut_convention)
    variants = resolve(args.variants)
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    suffix = "mockmatched" if args.reference_covariance == "mockmatched" else "data"
    payload: dict[str, object] = {}
    # The measured xi, loaded once: only the aliased variants and the guard
    # below read it, and both want the same two dicts gather() builds.
    gg = _load_xi(args.xi_reporting, args.nbins)
    gg_int = _load_xi(args.xi_integration, args.nbins_int)
    xi_kwargs = dict(min_sep=args.min_sep, max_sep=args.max_sep, nbins=args.nbins)

    with np.load(args.mock_modes, allow_pickle=False) as modes, tempfile.TemporaryDirectory() as tmp:
        staging = Path(tmp)
        # Each bundle is one variant; they unpack into a shared staging tree
        # because their member names are already variant-qualified.
        draws: dict[str, int] = {}
        for path in args.covariance_bundle:
            with tarfile.open(path, "r:gz") as bundle:
                bundle.extractall(staging)
            provenance = json.loads((staging / "provenance.json").read_text())
            for name in provenance["variants"]:
                draws[name] = int(provenance["n_samples"])
        available = {str(value) for value in np.asarray(modes["variants"])}
        for variant in variants:
            if variant.name not in available:
                raise ValueError(
                    f"{variant.name} is not in the mock pure-mode scan "
                    f"({sorted(available)}); the two stages must run the same "
                    "variant list."
                )
            # A variant with covariance_from set reads its donor's bundle:
            # its own propagated covariance is the same object, so building it
            # again would be a fake dependency. The alias is resolved here, in
            # the one place that maps a variant onto bundle members.
            source = variant.covariance_source
            reference = staging / f"{source}__{suffix}.npz"
            data_cov = staging / f"{source}__data.npz"
            for path in (reference, data_cov):
                if not path.is_file():
                    raise FileNotFoundError(
                        f"{path.name} missing from the covariance bundle"
                        + (
                            f" (needed for {variant.name} via covariance_from="
                            f"{variant.covariance_from!r})"
                            if variant.covariance_from
                            else ""
                        )
                    )
            if source not in draws:
                raise ValueError(
                    f"no covariance bundle supplied for {source}"
                    + (
                        f", the covariance source of {variant.name}"
                        if variant.covariance_from
                        else ""
                    )
                )
            # The Hartlap correction has to match the draws the covariance was
            # actually built from, not the pipeline's fiducial 2000.
            hartlap_n = draws[source]
            entry = _variant_statistics(modes, variant, reference, hartlap_n)

            if variant.covariance_from:
                # <source>__data.npz carries the SOURCE's data pure-mode
                # vectors as well as the covariance they share, and the data
                # vector is precisely what an aliased variant changes. Running
                # calculate_ptes on it as-is would report the donor's PTE under
                # this variant's name, so the vector is rebuilt here from the
                # measured xi -- through gather_pure_eb_chunks.data_vector, the
                # same function that produced the stored one -- and written
                # beside the borrowed covariance.
                #
                # GUARD: before trusting the rebuild, the donor's own vector is
                # recomputed the same way and checked against the one stored in
                # its bundle. Identical inputs through identical code must give
                # identical bytes; if they do not, the rebuild is not the
                # gather path and the PTE would be quietly wrong.
                donor = BY_NAME[variant.covariance_from]
                donor_modes = _rebuild_data_npz(
                    donor, data_cov, gg, gg_int, None, **xi_kwargs
                )
                with np.load(data_cov, allow_pickle=False) as stored:
                    mismatch = {
                        key: float(
                            np.max(np.abs(np.asarray(values, dtype=float) - ref))
                        )
                        for key, values in zip(_EB_KEYS, donor_modes)
                        for ref in [np.asarray(stored[key], dtype=float)]
                        # Not array_equal: the matvec's reduction order depends
                        # on the BLAS thread count, so a vector rebuilt on a
                        # different core count can differ in the last bits
                        # (measured at 4e-20 on values of order 1e-6, i.e. 1e-14
                        # relative). A wrong code path is not subtle -- the
                        # local-term substitution this guard protects moves the
                        # vector by O(1) relative -- so the tolerance is set
                        # eight orders below any real difference and two above
                        # the observed noise.
                        if not np.allclose(values, ref, rtol=1e-10, atol=1e-18)
                    }
                if mismatch:
                    raise ValueError(
                        f"in-stage data vector for {donor.name} does not "
                        f"reproduce its bundle's stored one (max|diff| per "
                        f"block: {mismatch}); {variant.name}'s rebuilt vector "
                        "cannot be trusted. The two must come from the same "
                        "gather code path."
                    )
                data_cov = staging / f"{variant.name}__data.npz"
                _rebuild_data_npz(
                    variant, staging / f"{source}__data.npz", gg, gg_int,
                    data_cov, **xi_kwargs
                )
                entry["data_pte_note"] = (
                    f"covariance borrowed from {variant.covariance_from} "
                    "(identical to round-off); the data vector is this "
                    "variant's own, rebuilt in-stage from the measured xi and "
                    f"checked against {donor.name}'s stored vector."
                )
            # The DATA's PTE under this variant, through the variant's own
            # covariance -- the same calculate_pure_eb_ptes call the paper's
            # single-variant output makes.
            pte_path = staging / f"{variant.name}__data_ptes.npz"
            calculate_ptes(
                version=f"scan_{variant.name}",
                blind="A",
                pure_eb_data=str(data_cov),
                cov_integration=args.cov_integration,
                npatch=args.npatch,
                out_file=str(pte_path),
            )
            with np.load(pte_path, allow_pickle=False) as ptes:
                entry["data_pte"] = {
                    name: _finite(ptes[key][window_min_bin, window_max_bin])
                    for name, key in DATA_PTE_KEYS.items()
                }
            entry.update(variant.as_dict())
            entry["reference_covariance"] = args.reference_covariance
            entry["covariance_draws"] = hartlap_n
            payload[variant.name] = entry
            print(f"{variant.name}: done", flush=True)

    theta_min, theta_max = PURE_THETA_BY_CONVENTION[args.scale_cut_convention]
    result = {
        "variant_order": [variant.name for variant in variants],
        "adopted": ADOPTED,
        "published": PUBLISHED,
        "dropped": DROPPED,
        "reference_covariance": args.reference_covariance,
        "covariance_draws": draws,
        "scale_cut_convention": args.scale_cut_convention,
        "theta_cut": [theta_min, theta_max],
        "window_bins": [window_min_bin, window_max_bin],
        "variants": payload,
    }
    with out_path.open("w") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Saved {out_path}", flush=True)

    header = f"{'':<26}" + "".join(f"{name:>14}" for name in result["variant_order"])
    print(header)
    rows = [
        ("mean-vector PTE xi+^B", lambda e: e["per_mock"]["xi+^B"]["mean_vector"]["pte"]),
        ("mean-vector PTE xi-^B", lambda e: e["per_mock"]["xi-^B"]["mean_vector"]["pte"]),
        ("KS p xi+^B", lambda e: e["per_mock"]["xi+^B"]["ks_pvalue"]),
        ("KS p xi-^B", lambda e: e["per_mock"]["xi-^B"]["ks_pvalue"]),
        ("KS p xi_tot^B", lambda e: e["per_mock"]["xi_tot^B"]["ks_pvalue"]),
        ("DATA PTE xi+^B", lambda e: e["data_pte"]["xi+^B"]),
        ("DATA PTE xi-^B", lambda e: e["data_pte"]["xi-^B"]),
        ("DATA PTE xi_tot^B", lambda e: e["data_pte"]["xi_tot^B"]),
    ]
    for label, getter in rows:
        cells = ""
        for name in result["variant_order"]:
            value = getter(payload[name])
            cells += f"{'--' if value is None else f'{value:.4g}':>14}"
        print(f"{label:<26}{cells}")


if __name__ == "__main__":
    main()
