"""Compute COSEBIS E_n/B_n for GLASS mock realizations.

`compute_one` holds the transform and the recipes call the batch CLI, whose
`--scale-cut-convention` carries the declared `scale_cut_convention` decision.
"""

import argparse
from pathlib import Path

import numpy as np
import treecorr
from cosmo_numba.B_modes.cosebis import COSEBIS

from sp_validation.b_modes import scale_cut_to_bins

import mock_scatter_pool

# The fiducial θ cut in its two forms.  `grid_bin_edges` is edges 9 and 16 of
# the reporting grid geomspace(1', 250', 21) — what the published pipeline cut
# at; `rounded_arcmin` is what the paper text prints.  Kept in step with
# `COSEBI_THETA_BY_CONVENTION` in mock_campaign_analysis.py and with
# papers/bmodes/scripts/scale_cuts.py, which resolves the same choice from the
# reporting_binning decision for the paper-facing scripts.
_REPORTING_EDGES = np.geomspace(1.0, 250.0, 21)
THETA_BY_CONVENTION = {
    "grid_bin_edges": (float(_REPORTING_EDGES[9]), float(_REPORTING_EDGES[16])),
    "rounded_arcmin": (12.0, 83.0),
}


class _Grid:
    """The three attributes `scale_cut_to_bins` reads off a GGCorrelation."""

    def __init__(self, meanr, left_edges, right_edges):
        self.meanr = meanr
        self.left_edges = left_edges
        self.right_edges = right_edges


def compute_from_grid(
    meanr, xip, xim, left_edges, right_edges, *, nmodes, theta_min, theta_max
):
    """COSEBI modes for one realisation's ξ± sampled on a fine grid."""

    grid = _Grid(meanr, left_edges, right_edges)
    start, stop = scale_cut_to_bins(grid, theta_min, theta_max)
    inds = np.arange(start, stop)
    # .astype(float) ensures native byte order (FITS big-endian → numba needs native)
    theta_cut = np.asarray(meanr)[inds].astype(float)
    xip_cut = np.asarray(xip)[inds].astype(float)
    xim_cut = np.asarray(xim)[inds].astype(float)

    cosebis_obj = COSEBIS(np.min(theta_cut), np.max(theta_cut), nmodes, precision=120)
    En, Bn = cosebis_obj.cosebis_from_xipm(theta_cut, xip_cut, xim_cut, parallel=True)
    return {
        "En": np.asarray(En, dtype=float),
        "Bn": np.asarray(Bn, dtype=float),
        "theta_min_actual": float(np.min(theta_cut)),
        "theta_max_actual": float(np.max(theta_cut)),
    }


def compute_one(xi_path, out_path, *, nmodes, theta_min, theta_max):
    """COSEBI modes for one mock's fine-binned ξ± FITS."""

    gg = treecorr.GGCorrelation(
        min_sep=0.5, max_sep=500, nbins=1000, sep_units="arcmin"
    )
    gg.read(xi_path)

    modes = compute_from_grid(
        gg.meanr,
        gg.xip,
        gg.xim,
        gg.left_edges,
        gg.right_edges,
        nmodes=nmodes,
        theta_min=theta_min,
        theta_max=theta_max,
    )
    np.savez(out_path, **modes)
    print(
        f"Saved {out_path}: {nmodes} modes, "
        f"θ=[{modes['theta_min_actual']:.2f}, {modes['theta_max_actual']:.2f}]"
    )
    return out_path


def _ensemble_worker(index):
    """One realisation, reading the ensemble the pool initializer shared."""

    state = mock_scatter_pool.state()
    return compute_from_grid(
        state["meanr"][index],
        state["xip"][index],
        state["xim"][index],
        state["left_edges"],
        state["right_edges"],
        nmodes=state["nmodes"],
        theta_min=state["theta_min"],
        theta_max=state["theta_max"],
    )


def compute_ensemble(ensemble_path, out_path, *, nmodes, theta_min, theta_max):
    """COSEBI modes for every realisation in a stacked ξ± ensemble."""

    with np.load(ensemble_path, allow_pickle=False) as data:
        mock_ids = np.asarray(data["mock_ids"])
        meanr, xip, xim = (np.asarray(data[key]) for key in ("meanr", "xip", "xim"))
        left_edges = np.asarray(data["left_edges"], dtype=float)
        right_edges = np.asarray(data["right_edges"], dtype=float)

    modes = mock_scatter_pool.map_realisations(
        _ensemble_worker,
        range(len(mock_ids)),
        {
            "meanr": meanr,
            "xip": xip,
            "xim": xim,
            "left_edges": left_edges,
            "right_edges": right_edges,
            "nmodes": nmodes,
            "theta_min": theta_min,
            "theta_max": theta_max,
        },
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        mock_ids=mock_ids,
        En=np.stack([entry["En"] for entry in modes]),
        Bn=np.stack([entry["Bn"] for entry in modes]),
        theta_min_actual=np.asarray(
            [entry["theta_min_actual"] for entry in modes], dtype=float
        ),
        theta_max_actual=np.asarray(
            [entry["theta_max_actual"] for entry in modes], dtype=float
        ),
        theta_min_requested=float(theta_min),
        theta_max_requested=float(theta_max),
    )
    print(
        f"Saved {out_path}: {len(mock_ids)} realisations x {nmodes} modes, "
        f"θ=[{modes[0]['theta_min_actual']:.2f}, "
        f"{modes[0]['theta_max_actual']:.2f}]",
        flush=True,
    )
    return out_path


def _mock_id_from_path(path):
    for token in Path(path).stem.split("_"):
        if token.isdigit() and len(token) == 5:
            return token
    raise ValueError(f"Cannot determine mock id from {path}")


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xi-ensemble",
        help="Stacked mock ξ± ensemble npz (mock_xi_ensemble.py); with --out-file.",
    )
    parser.add_argument(
        "--out-file",
        help="Exact npz path for the COSEBI mode ensemble.",
    )
    parser.add_argument(
        "--fine-dir",
        help="Directory of gg_glass_mock_<id>_nbins=1000.fits files (batch mode).",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for cosebis_glass_mock_<id>.npz (batch mode).",
    )
    parser.add_argument("--nmodes", type=int, default=20)
    parser.add_argument(
        "--scale-cut-convention",
        choices=tuple(THETA_BY_CONVENTION),
        help=(
            "How the fiducial COSEBI θ cut is read: 'grid_bin_edges' (edges 9 "
            "and 16 of the geomspace(1', 250', 21) reporting grid, what the "
            "published pipeline cut at) or 'rounded_arcmin' ([12, 83], what the "
            "paper text prints).  Mutually exclusive with the explicit "
            "--theta-min/--theta-max pair."
        ),
    )
    parser.add_argument("--theta-min", type=float)
    parser.add_argument("--theta-max", type=float)
    args = parser.parse_args(argv)
    explicit = args.theta_min is not None and args.theta_max is not None
    if bool(args.scale_cut_convention) == explicit:
        parser.error(
            "give either --scale-cut-convention or both of "
            "--theta-min/--theta-max, not both and not neither"
        )
    if args.scale_cut_convention:
        args.theta_min, args.theta_max = THETA_BY_CONVENTION[
            args.scale_cut_convention
        ]
    ensemble_mode = args.xi_ensemble is not None
    if ensemble_mode != (args.out_file is not None):
        parser.error("--xi-ensemble and --out-file must be given together")
    if not ensemble_mode and not (args.fine_dir and args.output_dir):
        parser.error("batch mode needs --fine-dir and --output-dir")
    return args


def _standalone_main(argv=None):
    args = _parser(argv)

    if args.xi_ensemble:
        compute_ensemble(
            args.xi_ensemble,
            args.out_file,
            nmodes=args.nmodes,
            theta_min=args.theta_min,
            theta_max=args.theta_max,
        )
        return

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fine_paths = sorted(Path(args.fine_dir).glob("gg_glass_mock_*_nbins=1000.fits"))
    if not fine_paths:
        raise FileNotFoundError(f"No fine ξ files under {args.fine_dir}")

    for fine_path in fine_paths:
        mock_id = _mock_id_from_path(fine_path)
        out_path = output_dir / f"cosebis_glass_mock_{mock_id}.npz"
        if out_path.is_file():
            continue
        compute_one(
            fine_path,
            out_path,
            nmodes=args.nmodes,
            theta_min=args.theta_min,
            theta_max=args.theta_max,
        )


def _snakemake_main(smk):
    compute_one(
        smk.input.xi,
        smk.output.cosebis,
        nmodes=smk.params.nmodes,
        theta_min=smk.params.theta_min,
        theta_max=smk.params.theta_max,
    )


if "snakemake" in globals():
    _snakemake_main(snakemake)  # noqa: F821
elif __name__ == "__main__":
    _standalone_main()
