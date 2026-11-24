"""Gather the per-row COSEBI PTE NPZs of the cosebis_pte_grid sweep into one NPZ.

`compute_cosebis_pte_single.py` loops the 206 stable (i_min, i_max) scale-cut
pairs and writes one gathered NPZ. That loop is the slowest step of the
cosebis analysis, so the `cosebis_pte_grid` astra-sweep splits it by i_min row
across nodes — each member runs the
same CLI with `--row i`, producing an NPZ of the *identical* shape and key order
whose `pairs` array holds only that row.

This is the gather half. It rebuilds the global pair list from
`_pte_scale_cut_pairs()` — the same function the members restrict — and places
each member's values by their (i_min, i_max) key, so the merged NPZ reproduces
the serial output exactly: same keys in the same order, same dtypes, same
shapes, same NaN placement. Nothing is recomputed and nothing is inferred from
NaN sentinels; a pair is copied from the one member that owns it, and coverage
is asserted to be exactly the global list, once each.

    python merge_cosebis_pte_rows.py --members <row.npz>... --out-file <path.npz>
"""

import argparse

import numpy as np

from compute_cosebis_pte_single import _pte_scale_cut_pairs

# The np.savez keyword order of compute_cosebis_pte_single.main(), split by how
# each group is merged. The concatenation of these three tuples, with "pairs"
# between the first and the second, IS the output contract.
SCALARS = ("version", "blind", "nmodes", "mode_subsets", "theta_grid")
PER_PAIR = ("theta_min", "theta_max", "En", "Bn", "convergence_failure")
MATRICES = (
    "pte_B_6",
    "chi2_B_6",
    "pte_E_6",
    "chi2_E_6",
    "pte_B_20",
    "chi2_B_20",
    "pte_E_20",
    "chi2_E_20",
)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--members",
        nargs="+",
        required=True,
        help="Per-row PTE NPZ paths (compute_cosebis_pte_single.py --row)",
    )
    ap.add_argument(
        "--out-file", required=True, help="Merged full-grid NPZ (lc {output})"
    )
    a = ap.parse_args(argv)

    members = [np.load(path) for path in a.members]
    ref = members[0]

    # The grid-wide metadata is version-independent and must agree exactly:
    # a disagreement means the members were not run against the same catalog,
    # theta grid or mode count, and merging them would be meaningless.
    for path, m in zip(a.members[1:], members[1:]):
        for key in SCALARS:
            if not np.array_equal(ref[key], m[key]):
                raise SystemExit(f"member {path} disagrees with {a.members[0]} on {key}")

    pairs = _pte_scale_cut_pairs()
    index = {(i, j): k for k, (i, j) in enumerate(pairs)}
    n_pairs = len(pairs)
    nmodes = int(ref["nmodes"])
    n_grid = int(ref["theta_grid"].shape[0])

    # Built in the serial np.savez keyword order; np.savez preserves **kwargs
    # order, so the merged archive carries the same entries in the same sequence.
    out = {key: ref[key] for key in SCALARS}
    out["pairs"] = np.array(pairs, dtype=int)
    out["theta_min"] = np.full(n_pairs, np.nan)
    out["theta_max"] = np.full(n_pairs, np.nan)
    out["En"] = np.full((n_pairs, nmodes), np.nan)
    out["Bn"] = np.full((n_pairs, nmodes), np.nan)
    out["convergence_failure"] = np.zeros(n_pairs, dtype=bool)
    for key in MATRICES:
        out[key] = np.full((n_grid, n_grid), np.nan)

    owner = {}
    for path, m in zip(a.members, members):
        for k_local, (i_min, i_max) in enumerate(m["pairs"]):
            key = (int(i_min), int(i_max))
            if key not in index:
                raise SystemExit(f"member {path} carries pair {key}, not in the grid")
            if key in owner:
                raise SystemExit(f"pair {key} produced by both {owner[key]} and {path}")
            owner[key] = path
            k = index[key]
            out["theta_min"][k] = m["theta_min"][k_local]
            out["theta_max"][k] = m["theta_max"][k_local]
            out["En"][k, :] = m["En"][k_local, :]
            out["Bn"][k, :] = m["Bn"][k_local, :]
            out["convergence_failure"][k] = m["convergence_failure"][k_local]
            for mat in MATRICES:
                out[mat][key] = m[mat][key]

    missing = [p for p in pairs if p not in owner]
    if missing:
        raise SystemExit(f"{len(missing)} pairs covered by no member: {missing}")

    with open(a.out_file, "wb") as f:
        np.savez(f, **out)
    print(
        f"Merged {len(a.members)} rows ({n_pairs} pairs, "
        f"{int(out['convergence_failure'].sum())} convergence failures) "
        f"into {a.out_file}"
    )


if __name__ == "__main__":
    main()
