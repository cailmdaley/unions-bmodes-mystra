"""Stack per-realisation TreeCorr ξ± files into one mock-ensemble array.

The gather of a `run_glass_mock_2pcf.py` sweep.  Each member is one mock's
GGCorrelation FITS; the ensemble is the realisation axis made explicit — a
single npz holding `(n_mocks, nbins)` blocks of `meanr`, `xip`, `xim`,
`npairs` and `weight` on a shared binning, plus the mock IDs.

This is what the campaign's downstream stages actually consume: none of them
wants one mock, all of them want the suite.  Keeping the ensemble as arrays
rather than a bundle of files means the consumers read numbers, not a
directory.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import treecorr

from decision_values import INTEGRATION_GRID

_ID_RE = re.compile(r"_(\d{5})(?:[_.]|$)")
_BLOCKS = ("meanr", "xip", "xim", "npairs", "weight")


def _mock_id_from_path(path: Path) -> str:
    match = _ID_RE.search(path.name)
    if match is None:
        raise ValueError(f"Could not recover a five-digit mock ID from {path}")
    return match.group(1)


def stack(members, *, min_sep, max_sep, nbins):
    """Read every member and return the ensemble arrays."""

    members = [Path(path) for path in members]
    if not members:
        raise ValueError("no ensemble members given")

    collected = {name: [] for name in _BLOCKS}
    mock_ids = []
    edges = None
    for path in sorted(members, key=_mock_id_from_path):
        gg = treecorr.GGCorrelation(
            min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units="arcmin"
        )
        gg.read(str(path))
        if len(gg.meanr) != nbins:
            raise ValueError(
                f"{path} has {len(gg.meanr)} bins, expected {nbins}"
            )
        member_edges = (
            np.asarray(gg.left_edges, dtype=float),
            np.asarray(gg.right_edges, dtype=float),
        )
        if edges is None:
            edges = member_edges
        elif not all(
            np.allclose(a, b, rtol=1e-12) for a, b in zip(edges, member_edges)
        ):
            raise ValueError(f"{path} does not share the ensemble binning")
        for name in _BLOCKS:
            collected[name].append(np.asarray(getattr(gg, name), dtype=float))
        mock_ids.append(_mock_id_from_path(path))

    return {
        "mock_ids": np.asarray(mock_ids),
        "left_edges": edges[0],
        "right_edges": edges[1],
        "min_sep": float(min_sep),
        "max_sep": float(max_sep),
        "nbins": int(nbins),
        **{name: np.stack(values) for name, values in collected.items()},
    }


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--members", nargs="+", required=True)
    parser.add_argument(
        "--integration-grid", type=INTEGRATION_GRID.__getitem__, required=True
    )
    parser.add_argument("--out-file", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = _parser(argv)
    ensemble = stack(
        args.members,
        min_sep=args.integration_grid.min_sep,
        max_sep=args.integration_grid.max_sep,
        nbins=args.integration_grid.nbins,
    )
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, **ensemble)
    print(
        f"Saved {out_path}: {len(ensemble['mock_ids'])} realisations "
        f"x {args.integration_grid.nbins} bins",
        flush=True,
    )


if __name__ == "__main__":
    main()
