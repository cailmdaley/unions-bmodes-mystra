"""Compute the pure E/B configuration-space data vector for GLASS mocks.

The reporting grid is the external 20-bin GLASS ``xi`` product (1--250
arcmin).  The integration grid is the locally-produced 1000-bin ``gg``
product (0.5--300 arcmin).  This is deliberately only the data-vector part of
``sp_validation.b_modes.calculate_pure_eb_correlation``: no covariance or
jackknife machinery is run here.  The transform's quadrature is selected by
``--integration-method`` (the ``pure_eb_integration_method`` decision), which
also selects the geometry: ``adaptive_dqags`` transforms pointwise on the
reporting grid, ``fixed_operator`` transforms on the fine grid and averages the
pure modes into the reporting bins under the weights ``--rebinning-weights``
(the ``pure_eb_rebinning_weights`` decision) selects.

The CLI accepts either explicit mock IDs or a glob of fine-``xi`` files.
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np
import treecorr

import mock_scatter_pool

# The two script trees are separate sys.path roots; the transform helper that
# resolves `pure_eb_integration_method` lives with the paper-facing scripts and
# is shared rather than duplicated.
sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "papers" / "bmodes" / "scripts")
)
from pure_eb_transform import (  # noqa: E402
    AveragedPureEBTransform,
    PureEBTransform,
    add_integration_limits_argument,
    add_integration_method_argument,
    add_rebinning_weights_argument,
    check_rebinning_weights,
    geometry_for,
    resolve_integration_limits,
)

_EB_KEYS = ("xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")
_ID_RE = re.compile(r"(?:glass_mock_|mock_)(\d{5})")


def _native(values):
    """Return a native-byte-order float array for numba inputs."""

    return np.asarray(values, dtype=float)


def _print_drop_counts(counts, mock_id):
    details = ", ".join(f"{name}={count}" for name, count in zip(_EB_KEYS, counts))
    print(
        f"Dropped non-finite outputs for mock {mock_id}: {details}",
        flush=True,
    )


def _mock_id_from_path(path: Path) -> str:
    match = _ID_RE.search(path.name)
    if match is None:
        raise ValueError(f"Could not recover a five-digit mock ID from {path}")
    return match.group(1)


def _reporting_path(reporting_dir: Path, mock_id: str) -> Path:
    """External GLASS products use `_4096_nbins=20`; locally produced
    reporting ξ (v2/masked) use `_nbins=20`.  Prefer the external
    name, fall back to the local one."""
    external = reporting_dir / f"xi_glass_mock_{mock_id}_4096_nbins=20.fits"
    if external.is_file():
        return external
    return reporting_dir / f"xi_glass_mock_{mock_id}_nbins=20.fits"


def _output_path(output_dir: Path, mock_id: str) -> Path:
    return output_dir / f"pure_eb_glass_mock_{mock_id}.npz"


class _Grid:
    """The reporting-grid attributes the transform and binning matrix read."""

    def __init__(self, meanr, xip, xim, npairs, left_edges, right_edges):
        self.meanr = meanr
        self.xip = xip
        self.xim = xim
        self.npairs = npairs
        self.left_edges = left_edges
        self.right_edges = right_edges


def compute_from_grids(
    mock_id: str,
    gg_int: _Grid,
    gg: _Grid,
    *,
    integration_method: str,
    rebinning_weights: str,
    integration_limits: str,
    transform_averaged: bool = False,
    geometry: str | None = None,
    local_from_eval_grid: bool = False,
) -> dict:
    """Transform one realisation's fine/reporting ξ± pair into pure modes."""

    # Match calculate_pure_eb_correlation: tmin/tmax are the actual reporting
    # bin edges, not the nominal constructor arguments.  Explicit casts also
    # handle FITS big-endian arrays before they enter numba.
    tmin = float(gg.left_edges[0])
    tmax = float(gg.right_edges[-1])
    theta_int = _native(gg_int.meanr)
    xip_int = _native(gg_int.xip)
    xim_int = _native(gg_int.xim)
    npairs_int = _native(gg_int.npairs)
    # transform_averaged forces the averaged geometry; `geometry` is the explicit
    # override. Resolved through geometry_for so the method's default still
    # applies when neither is given -- an OR here would make fixed_operator
    # mean "averaged" unconditionally and put pointwise_fixedop out of reach.
    geometry = geometry_for(
        integration_method, "averaged" if transform_averaged else geometry
    )
    averaged = geometry == "averaged"
    theta_window = resolve_integration_limits(integration_limits, geometry)
    check_rebinning_weights(rebinning_weights, geometry)
    if averaged:
        transform = AveragedPureEBTransform(
            integration_method,
            left_edges=_native(gg.left_edges),
            right_edges=_native(gg.right_edges),
            theta_int=theta_int,
            rebinning_weights=rebinning_weights,
            npairs=npairs_int,
            theta_window=theta_window,
            parallel=True,
            # The averaged hybrid: same operator and same uniform average for
            # the integral term, but the local 0.5*(xi_p +/- xi_m) term taken
            # from the 20-bin reporting measurement instead of the fine grid's
            # own bin average of it.
            local_from_eval_grid=local_from_eval_grid,
        )
        pure_modes = transform(
            xip_int, xim_int, _native(gg.xip), _native(gg.xim)
        )
        _print_drop_counts(transform.last_dropped, mock_id)
        snap_offset_max = 0.0
    else:
        if theta_window is not None:
            lo, hi = theta_window
            keep = np.ones(theta_int.size, dtype=bool)
            if lo is not None:
                keep &= theta_int >= float(lo)
            if hi is not None:
                keep &= theta_int <= float(hi)
            if np.count_nonzero(keep) < 10:
                raise ValueError(
                    f"theta_window {theta_window} leaves too few fine xi bins: "
                    f"{np.count_nonzero(keep)}"
                )
            theta_int, xip_int, xim_int, npairs_int = (
                values[keep] for values in (theta_int, xip_int, xim_int, npairs_int)
            )
        transform = PureEBTransform(
            integration_method,
            theta=_native(gg.meanr),
            theta_int=theta_int,
            tmin=tmin,
            tmax=tmax,
            parallel=True,
            # The single-grid operator cannot evaluate between nodes, so the
            # pointwise operator chain evaluates at the nearest fine node in
            # log theta instead. dqags interpolates and needs no snap, so this
            # leaves the published pointwise chain untouched.
            snap_eval=integration_method == "fixed_operator",
            # The hybrid corner: snap for the integral term, but read the local
            # 0.5*(xi_p +/- xi_m) term from the reporting measurement rather
            # than from the single snapped fine node. A no-op unless snapping
            # ran, and refused under adaptive_dqags, which already does this.
            local_from_eval_grid=local_from_eval_grid,
        )
        pure_modes = transform(_native(gg.xip), _native(gg.xim), xip_int, xim_int)
        snap_offset_max = transform.snap_offset_max
    arrays = dict(zip(_EB_KEYS, (_native(values) for values in pure_modes)))

    if not all(np.all(np.isfinite(values)) for values in arrays.values()):
        raise ValueError(f"Non-finite pure E/B output for mock {mock_id}")

    return {
        "mock_id": np.asarray(mock_id),
        "theta": _native(gg.meanr),
        "reporting_left_edges": _native(gg.left_edges),
        "reporting_right_edges": _native(gg.right_edges),
        "tmin": tmin,
        "tmax": tmax,
        "integration_limits": integration_limits,
        "transform": "averaged" if averaged else "pointwise",
        "rebinning_weights": rebinning_weights,
        "local_from_eval_grid": bool(local_from_eval_grid),
        # Largest |theta_snapped/theta_requested - 1| over the evaluation
        # points; 0.0 for every chain that evaluates where it was asked to.
        "snap_offset_max": snap_offset_max,
        **arrays,
    }


def _grid_from_fits(path, *, min_sep, max_sep, nbins) -> _Grid:
    gg = treecorr.GGCorrelation(
        min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units="arcmin"
    )
    gg.read(str(path))
    return _Grid(gg.meanr, gg.xip, gg.xim, gg.npairs, gg.left_edges, gg.right_edges)


def _grid_from_ensemble(data, index) -> _Grid:
    return _Grid(
        data["meanr"][index],
        data["xip"][index],
        data["xim"][index],
        data["npairs"][index],
        data["left_edges"],
        data["right_edges"],
    )


def compute_one(
    mock_id: str,
    fine_xi_path: str | Path,
    reporting_xi_path: str | Path,
    output_path: str | Path,
    *,
    integration_method: str,
    rebinning_weights: str,
    integration_limits: str,
    geometry: str | None = None,
    fine_min_sep: float = 0.5,
    fine_max_sep: float = 300.0,
    fine_nbins: int = 1000,
    reporting_min_sep: float = 1.0,
    reporting_max_sep: float = 250.0,
    reporting_nbins: int = 20,
    transform_averaged: bool = False,
    local_from_eval_grid: bool = False,
) -> Path:
    """Transform one fine/reporting FITS pair and write its pure-mode arrays."""

    fine_xi_path = Path(fine_xi_path)
    reporting_xi_path = Path(reporting_xi_path)
    output_path = Path(output_path)
    if not fine_xi_path.is_file():
        raise FileNotFoundError(fine_xi_path)
    if not reporting_xi_path.is_file():
        raise FileNotFoundError(reporting_xi_path)

    result = compute_from_grids(
        mock_id,
        _grid_from_fits(
            fine_xi_path,
            min_sep=fine_min_sep,
            max_sep=fine_max_sep,
            nbins=fine_nbins,
        ),
        _grid_from_fits(
            reporting_xi_path,
            min_sep=reporting_min_sep,
            max_sep=reporting_max_sep,
            nbins=reporting_nbins,
        ),
        integration_method=integration_method,
        rebinning_weights=rebinning_weights,
        integration_limits=integration_limits,
        transform_averaged=transform_averaged,
        geometry=geometry,
        local_from_eval_grid=local_from_eval_grid,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output_path, **result)
    print(
        f"Saved {output_path}: mock={mock_id}, ntheta={len(result['theta'])}, "
        f"theta=[{result['tmin']:.6g}, {result['tmax']:.6g}] arcmin",
        flush=True,
    )
    return output_path


def _ensemble_worker(index):
    """One realisation, reading the ensembles the pool initializer shared."""

    state = mock_scatter_pool.state()
    return compute_from_grids(
        state["mock_ids"][index],
        _grid_from_ensemble(state["fine"], index),
        _grid_from_ensemble(state["reporting"], index),
        integration_method=state["integration_method"],
        rebinning_weights=state["rebinning_weights"],
        integration_limits=state["integration_limits"],
        transform_averaged=state["transform_averaged"],
        geometry=state["geometry"],
        local_from_eval_grid=state["local_from_eval_grid"],
    )


def compute_ensemble(
    fine_ensemble_path,
    reporting_ensemble_path,
    output_path,
    *,
    integration_method: str,
    rebinning_weights: str,
    integration_limits: str,
    transform_averaged: bool = False,
    geometry: str | None = None,
    local_from_eval_grid: bool = False,
) -> Path:
    """Transform every realisation of a fine/reporting ensemble pair."""

    with np.load(fine_ensemble_path, allow_pickle=False) as fine, np.load(
        reporting_ensemble_path, allow_pickle=False
    ) as reporting:
        fine_ids = [str(value) for value in np.asarray(fine["mock_ids"])]
        reporting_ids = [str(value) for value in np.asarray(reporting["mock_ids"])]
        # The reporting ensemble stacks the whole suite; the fine ensemble is
        # the staged subset actually measured. Align reporting to fine by ID.
        missing = sorted(set(fine_ids) - set(reporting_ids))
        if missing:
            raise ValueError(
                "reporting ensemble lacks realisations present in the fine "
                f"ensemble: {missing[:10]}"
            )
        reporting_index = np.asarray(
            [reporting_ids.index(mock_id) for mock_id in fine_ids]
        )
        # Materialise before forking: the workers share these arrays by fork
        # rather than pickling a lazy npz handle per realisation.
        keys = ("meanr", "xip", "xim", "npairs", "left_edges", "right_edges")
        per_realisation = ("meanr", "xip", "xim", "npairs")
        state = {
            "mock_ids": fine_ids,
            "fine": {key: np.asarray(fine[key]) for key in keys},
            "reporting": {
                key: np.asarray(reporting[key])[reporting_index]
                if key in per_realisation
                else np.asarray(reporting[key])
                for key in keys
            },
            "integration_method": integration_method,
            "rebinning_weights": rebinning_weights,
            "integration_limits": integration_limits,
            "transform_averaged": transform_averaged,
            "geometry": geometry,
            "local_from_eval_grid": local_from_eval_grid,
        }
    results = mock_scatter_pool.map_realisations(
        _ensemble_worker, range(len(fine_ids)), state
    )

    first = results[0]
    stacked = {key: np.stack([r[key] for r in results]) for key in _EB_KEYS}
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        mock_ids=np.asarray(fine_ids),
        theta=np.stack([r["theta"] for r in results]),
        reporting_left_edges=first["reporting_left_edges"],
        reporting_right_edges=first["reporting_right_edges"],
        tmin=first["tmin"],
        tmax=first["tmax"],
        integration_limits=first["integration_limits"],
        transform=first["transform"],
        rebinning_weights=first["rebinning_weights"],
        local_from_eval_grid=first["local_from_eval_grid"],
        snap_offset_max=first["snap_offset_max"],
        **stacked,
    )
    print(
        f"Saved {output_path}: {len(fine_ids)} realisations, "
        f"ntheta={first['theta'].size}, "
        f"theta=[{first['tmin']:.6g}, {first['tmax']:.6g}] arcmin",
        flush=True,
    )
    return output_path


def _standalone_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fine-ensemble",
        help="Stacked fine ξ± ensemble npz; with --reporting-ensemble/--out-file.",
    )
    parser.add_argument(
        "--reporting-ensemble",
        help="Stacked 20-bin reporting ξ± ensemble npz.",
    )
    parser.add_argument(
        "--out-file",
        help="Exact npz path for the pure-mode ensemble.",
    )
    parser.add_argument(
        "--mock-ids",
        nargs="+",
        help="Five-digit mock IDs; if omitted, IDs are recovered from --fine-glob.",
    )
    parser.add_argument(
        "--fine-glob",
        default="results/glass_mock/gg_glass_mock_*_nbins=1000.fits",
        help="Glob for fine mock ξ files (quoted so the script expands it).",
    )
    parser.add_argument(
        "--reporting-dir",
        default="/n09data/guerrini/glass_mock_v1.4.6/results",
        help="Directory containing external 20-bin GLASS ξ FITS files.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/glass_mock",
        help="Directory for pure_eb_glass_mock_<id>.npz files.",
    )
    add_integration_limits_argument(parser)
    parser.add_argument(
        "--geometry",
        choices=("pointwise", "averaged"),
        default=None,
        help=(
            "Override the geometry the integration method defaults to. "
            "'pointwise' with fixed_operator snaps the reporting abscissae "
            "onto the fine grid; see pure_eb_transform.geometry_for."
        ),
    )
    parser.add_argument(
        "--transform-averaged",
        action="store_true",
        help=(
            "Force the averaged geometry under adaptive_dqags; fixed_operator "
            "already implies it."
        ),
    )
    parser.add_argument(
        "--local-from-eval-grid",
        action="store_true",
        help=(
            "Pointwise fixed_operator only: take the local 0.5*(xi_p +/- xi_m) "
            "term from the reporting-grid xi at the requested abscissa rather "
            "than from the snapped fine node, as the dqags pointwise chain "
            "does. The integral term still comes from the fine grid at the "
            "snapped node."
        ),
    )
    add_integration_method_argument(parser)
    add_rebinning_weights_argument(parser)
    return parser.parse_args(argv)


def _standalone_main(argv=None):
    args = _standalone_args(argv)

    ensemble_flags = (args.fine_ensemble, args.reporting_ensemble, args.out_file)
    if any(ensemble_flags):
        if not all(ensemble_flags):
            raise SystemExit(
                "ensemble mode needs --fine-ensemble, --reporting-ensemble "
                "and --out-file together"
            )
        compute_ensemble(
            args.fine_ensemble,
            args.reporting_ensemble,
            args.out_file,
            integration_method=args.integration_method,
            rebinning_weights=args.rebinning_weights,
            integration_limits=args.integration_limits,
            transform_averaged=args.transform_averaged,
            geometry=args.geometry,
            local_from_eval_grid=args.local_from_eval_grid,
        )
        return

    fine_paths = {}
    if args.mock_ids:
        for raw_id in args.mock_ids:
            mock_id = str(raw_id).zfill(5)
            fine_path = Path(args.fine_glob.format(mock_id=mock_id))
            # An explicit ID may use a literal path pattern, or the default
            # wildcard.  Resolve both forms without importing glob globally.
            if any(char in args.fine_glob for char in "*?["):
                matches = sorted(Path(path) for path in glob.glob(args.fine_glob))
                matches = [p for p in matches if _mock_id_from_path(p) == mock_id]
                if matches:
                    fine_path = matches[0]
            if fine_path.is_file():
                fine_paths[mock_id] = fine_path
            else:
                print(f"Skipping {mock_id}: missing fine ξ {fine_path}")
    else:
        fine_paths = {
            _mock_id_from_path(path): path
            for path in sorted(Path(path) for path in glob.glob(args.fine_glob))
        }

    if not fine_paths:
        raise FileNotFoundError(f"No fine ξ files matched {args.fine_glob!r}")

    for mock_id, fine_path in sorted(fine_paths.items()):
        reporting_path = _reporting_path(Path(args.reporting_dir), mock_id)
        output_path = _output_path(Path(args.output_dir), mock_id)
        # A concurrently-produced mock is allowed to be absent on the
        # external reporting side; leave it for a later invocation.
        if not reporting_path.is_file():
            print(f"Skipping {mock_id}: missing reporting ξ {reporting_path}")
            continue
        if output_path.is_file():
            continue
        compute_one(
            mock_id,
            fine_path,
            reporting_path,
            output_path,
            integration_method=args.integration_method,
            rebinning_weights=args.rebinning_weights,
            integration_limits=args.integration_limits,
            transform_averaged=args.transform_averaged,
            geometry=args.geometry,
            local_from_eval_grid=args.local_from_eval_grid,
        )


def _snakemake_main(smk):
    # The published chain's settings; recipes select the decision through the CLI.
    compute_one(
        str(smk.wildcards.mock_id),
        smk.input.fine,
        smk.input.reporting,
        smk.output.pure_eb,
        integration_method="adaptive_dqags",
        rebinning_weights="none",
        integration_limits="reporting_extent",
    )


# The __main__ guard on the CLI branch is what lets the scan stages import
# compute_ensemble instead of reimplementing it -- without it, importing this
# module parses argv and runs a job.
if "snakemake" in globals():
    _snakemake_main(snakemake)
elif __name__ == "__main__":
    _standalone_main()
