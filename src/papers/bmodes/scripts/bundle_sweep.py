"""Bundle astra-sweep member artifacts into the July-canonical sweep tar.gz.

The Design-B version sweeps (ξ±, pseudo-Cℓ, CosmoCov) are expressed as
astra-sweep generated sub-analyses: one lc output per (version, …) member, plus
a single *gather* output that is the public surface the hand-authored spec
references. Downstream consumers (``cosebis_version_comparison --results-dir``,
``cl_version_comparison --results-dir``, the pure E/B sweep ``--xi-sweep-dir``,
``--cov-dir``, ``config_space_pte_matrices --cosebis-pte-dir``) still read the
sweep as a *directory of canonically-named files*,
so the gather must reassemble exactly the July filenames.

lc names a member artifact ``results/<universe>/<output_id>.<format>``, which
carries no version/grid information, so this script re-derives the
member-id → canonical-filename map from the same matrix the ``.sweep.yaml``
declares (version list from :mod:`sweep_versions`, grids/bins/products
hard-coded per layout below) and applies the astra-sweep identifier munging
rule to match ``--names`` against it. Lookup is by name, not by position, and
an unmatched member id or a missing expected member is a hard error — so a
matrix change that this file does not track fails loudly instead of silently
bundling the wrong thing.

    python bundle_sweep.py --layout xi \
        --config .../config.yaml \
        --names <member_id>... --members <member_path>... \
        --out-file <path.tar.gz>
"""

import argparse
import os
import re
import shutil
import tarfile
import tempfile

import yaml

from sweep_versions import nonfiducial_versions
from paper_config import load_config

BLIND = "A"
BINNING = "powspace"

# ξ± grids (run_xi_sweep.GRIDS), npatch fixed at the paper fiducial.
XI_GRIDS = {
    "reporting": dict(min_sep=1.0, max_sep=250.0, nbins=20),
    "integration": dict(min_sep=0.5, max_sep=300.0, nbins=1000),
}
XI_NPATCH = 1

# pseudo-Cℓ: 32-bin reporting grid for every version; 96-bin fine grid only for
# the leak-corrected versions (it exists to feed the harmonic COSEBI scan).
CL_NBINS = (32, 96)
CL_PRODUCTS = ("cl", "cl_cov")

# CosmoCov integration grid (Gaussian-only, masked, blind A).
COV_GRID = dict(min_sep=0.5, max_sep=300.0, nbins=1000)

# The COSEBI PTE-matrix scan has no grid axis: one gathered NPZ per version.


def munge_id(value):
    """astra-sweep's identifier munging: lowercase, '.'→'p', other runs → '_'."""
    s = str(value).strip().lower()
    if s.startswith("-"):
        s = "neg" + s[1:]
    s = s.replace(".", "p")
    s = re.sub(r"[^a-z0-9_]+", "_", s)
    return s.strip("_")


def _xi_map(versions):
    out = {}
    for ver in versions:
        for grid, g in XI_GRIDS.items():
            member = f"xi_{munge_id(ver)}_{grid}"
            out[member] = (
                f"{ver}_xi_minsep={g['min_sep']}_maxsep={g['max_sep']}"
                f"_nbins={g['nbins']}_npatch={XI_NPATCH}.txt"
            )
    return out


def _cl_map(versions):
    out = {}
    for ver in versions:
        for nbins in CL_NBINS:
            if nbins == 96 and "_leak_corr" not in ver:
                continue  # fine grid is leak-corrected-only (harmonic COSEBI scan)
            for product in CL_PRODUCTS:
                member = f"pseudo_{product}_{munge_id(ver)}_nbins{nbins}"
                out[member] = (
                    f"pseudo_{product}_{ver}_blind={BLIND}_{BINNING}_nbins={nbins}.fits"
                )
    return out


def _cov_map(versions):
    out = {}
    for ver in versions:
        base = (
            f"covariance_{ver}_{BLIND}_g_minsep={COV_GRID['min_sep']}"
            f"_maxsep={COV_GRID['max_sep']}_nbins={COV_GRID['nbins']}_masked"
        )
        out[f"cov_{munge_id(ver)}"] = os.path.join(base, f"{base}_processed.txt")
    return out


def _cosebis_ptes_map(versions):
    out = {}
    for ver in versions:
        out[f"cosebis_ptes_{munge_id(ver)}"] = f"cosebis_ptes_{ver}_{BLIND}.npz"
    return out


LAYOUTS = {
    "xi": _xi_map,
    "cl": _cl_map,
    "cov": _cov_map,
    "cosebis_ptes": _cosebis_ptes_map,
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--layout", required=True, choices=sorted(LAYOUTS))
    ap.add_argument("--config", required=True, help="Absolute path to bmodes config.yaml")
    ap.add_argument("--names", nargs="+", required=True, help="Member output ids")
    ap.add_argument("--members", nargs="+", required=True, help="Member artifact paths")
    ap.add_argument("--out-file", required=True, help="Bundle path (lc {output})")
    a = ap.parse_args(argv)

    if len(a.names) != len(a.members):
        ap.error(
            f"--names ({len(a.names)}) and --members ({len(a.members)}) must be "
            "the same length; they come from the same expansion"
        )

    config = load_config(a.config)
    expected = LAYOUTS[a.layout](nonfiducial_versions(config))

    got = dict(zip(a.names, a.members))
    unknown = sorted(set(got) - set(expected))
    missing = sorted(set(expected) - set(got))
    if unknown or missing:
        ap.error(
            f"layout {a.layout!r} does not match the sweep matrix — "
            f"unknown members: {unknown}; missing members: {missing}. "
            "Update the layout table in bundle_sweep.py alongside the "
            ".sweep.yaml matrix."
        )

    scratch = tempfile.mkdtemp(prefix="bundle_sweep_", dir=os.environ.get("TMPDIR"))
    try:
        for member, canonical in sorted(expected.items()):
            dest = os.path.join(scratch, canonical)
            os.makedirs(os.path.dirname(dest) or scratch, exist_ok=True)
            shutil.copyfile(got[member], dest)
            print(f"[bundle_sweep] {member} -> {canonical}")
        with tarfile.open(a.out_file, "w:gz") as tar:
            tar.add(scratch, arcname=".")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    print(f"[bundle_sweep] wrote {a.out_file} ({len(expected)} entries)")


if __name__ == "__main__":
    main()
