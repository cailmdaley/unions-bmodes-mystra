"""Re-measure one GLASS realisation and difference it against the ingested one.

The fine ξ_± ensemble is read from the 21–23 August 2026 snakemake measurement
rather than remade — 5 000 core-hours the campaign does not have. This is the
receipt on that seam: one realisation measured here, under this spec's
environment and this spec's TreeCorr call, against the stored vector standing
in for it in the ensemble.

Writes a JSON verdict. A mismatch means the ingest is no longer standing in for
what this spec would compute, and the ingested suite is what to suspect.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import treecorr

# The measurement itself must be the one the sweep would have run, not a copy:
# import it from the script beside this one.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_glass_mock_2pcf import compute_one  # noqa: E402


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-dir", required=True)
    parser.add_argument("--xi-ensemble", required=True)
    parser.add_argument("--mock-id", default="00001")
    parser.add_argument("--min-sep", type=float, required=True)
    parser.add_argument("--max-sep", type=float, required=True)
    parser.add_argument("--nbins", type=int, required=True)
    parser.add_argument("--num-threads", type=int, default=None)
    parser.add_argument("--out-file", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = _parser(argv)
    mock_id = str(args.mock_id).zfill(5)

    ensemble = np.load(args.xi_ensemble)
    ids = [str(m) for m in ensemble["mock_ids"]]
    if mock_id not in ids:
        raise SystemExit(f"{mock_id} is not in {args.xi_ensemble}")
    row = ids.index(mock_id)

    num_threads = args.num_threads
    if num_threads is None:
        num_threads = int(os.environ.get("OMP_NUM_THREADS", "1") or 1)

    with tempfile.TemporaryDirectory() as tmp:
        fresh_path = Path(tmp) / f"gg_{mock_id}.fits"
        compute_one(
            Path(args.catalog_dir) / f"unions_glass_sim_{mock_id}_4096.fits",
            fresh_path,
            min_sep=args.min_sep,
            max_sep=args.max_sep,
            nbins=args.nbins,
            num_threads=num_threads,
        )
        gg = treecorr.GGCorrelation(
            min_sep=args.min_sep,
            max_sep=args.max_sep,
            nbins=args.nbins,
            sep_units="arcmin",
        )
        gg.read(str(fresh_path))

    report = {
        "mock_id": mock_id,
        "nbins": int(args.nbins),
        "num_threads": int(num_threads),
        "ingested_ensemble": str(args.xi_ensemble),
    }

    for name in ("meanr", "xip", "xim", "npairs", "weight"):
        fresh = np.asarray(getattr(gg, name), dtype=float)
        stored = np.asarray(ensemble[name][row], dtype=float)
        diff = fresh - stored
        scale = np.max(np.abs(stored))
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(stored != 0, np.abs(diff) / np.abs(stored), 0.0)
        report[name] = {
            "max_abs_diff": float(np.max(np.abs(diff))),
            "max_abs_diff_over_max_stored": float(np.max(np.abs(diff)) / scale)
            if scale
            else None,
            "max_fractional_diff": float(np.max(frac)),
            "identical": bool(np.array_equal(fresh, stored)),
        }

    edges_match = bool(
        np.allclose(ensemble["left_edges"], gg.left_edges, rtol=1e-12)
        and np.allclose(ensemble["right_edges"], gg.right_edges, rtol=1e-12)
    )
    report["bin_edges_match"] = edges_match
    report["verdict"] = (
        "ingest reproduces this spec's measurement"
        if edges_match
        and report["xip"]["max_fractional_diff"] < 1e-6
        and report["xim"]["max_fractional_diff"] < 1e-6
        else "MISMATCH - the ingested suite does not reproduce"
    )

    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
