"""Bundle astra-sweep member artifacts into one sweep tar.gz of canonical filenames.

The version sweeps (ξ±, Cℓ, CosmoCov) are astra-sweep generated
sub-analyses: one output per (version, …) member, plus a single *gather* output
that the hand-authored spec references. Downstream consumers (``cosebis_version_comparison --results-dir``,
``cl_version_comparison --results-dir``, the pure E/B sweep ``--xi-sweep-dir``,
``--cov-dir``, ``config_space_pte_matrices --cosebis-pte-dir``) read the
sweep as a *directory of canonically-named files*, so the gather reassembles
those filenames.

A member artifact is named ``results/<universe>/<output_id>.<format>``, which
carries no version/grid information, so this script re-derives the
member-id → canonical-filename map from the same matrix the ``.sweep.yaml``
declares (version list from :mod:`catalog_versions`, grids/bins/products
hard-coded per layout below) and applies the astra-sweep identifier munging
rule to match ``--names`` against it. Lookup is by name, not by position, and
an unmatched member id or a missing expected member is a hard error — so a
matrix change that this file does not track fails loudly instead of silently
bundling the wrong thing.

    python bundle_sweep.py --layout xi \
        --fiducial-catalog size_cut --leakage-correction applied \
        --names <member_id>... --members <member_path>... \
        --out-file <path.tar.gz>
"""

import argparse
import os
import re
import shutil
import tarfile
import tempfile

from catalog_versions import add_catalog_arguments, catalogs_from
from decision_values import HARMONIC_VALIDATION_BINNING

BLIND = "A"
BINNING = "powspace"

# ξ± grids, npatch fixed at the paper fiducial.
XI_GRIDS = {
    "reporting": dict(min_sep=1.0, max_sep=250.0, nbins=20),
    "integration": dict(min_sep=0.5, max_sep=300.0, nbins=1000),
}
XI_NPATCH = 1

# Cℓ: 32-bin reporting grid for every version; 96-bin fine grid only for
# the leak-corrected versions (it exists to feed the harmonic COSEBI scan).
CL_GRIDS = ("nbins32", "fine")
CL_PRODUCTS = ("cl", "cl_cov")

# CosmoCov integration grid (selected class, Gaussian area-only, blind A).
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


def _slot(ver):
    return "fiducial_counterpart" if ver == "SP_v1.4.6.3" else munge_id(ver)


def _xi_map(versions, catalogs):
    out = {}
    for ver in versions:
        for grid, g in XI_GRIDS.items():
            member = f"xi_{_slot(ver)}_{grid}"
            actual = catalogs.scan_version(ver)
            out[member] = (
                f"{actual}_xi_minsep={g['min_sep']}_maxsep={g['max_sep']}"
                f"_nbins={g['nbins']}_npatch={XI_NPATCH}.txt"
            )
    return out


def _cl_map(versions, fine_nbins, catalogs):
    out = {}
    for ver in versions:
        for grid in CL_GRIDS:
            fine = grid == "fine"
            if fine and "_leak_corr" not in ver:
                continue  # fine grid is leak-corrected-only (harmonic COSEBI scan)
            for product in CL_PRODUCTS:
                slot = munge_id(ver.removesuffix("_leak_corr")) if fine else _slot(ver)
                member = f"{product}_{slot}_{grid}"
                actual = catalogs.scan_version(ver, fine=fine)
                nbins = fine_nbins if fine else 32
                out[member] = (
                    f"{product}_{actual}_blind={BLIND}_{BINNING}_nbins={nbins}.fits"
                )
    return out


def _cov_map(versions, catalogs):
    out = {}
    for ver in versions:
        actual = catalogs.scan_version(ver)
        base = f"covariance_{actual}_{BLIND}_integration"
        out[f"cov_{_slot(ver)}"] = os.path.join(base, f"{base}_processed.txt")
    return out


def _cosebis_ptes_map(versions, catalogs):
    out = {}
    for ver in versions:
        actual = catalogs.scan_version(ver)
        out[f"cosebis_ptes_{_slot(ver)}"] = f"cosebis_ptes_{actual}_{BLIND}.npz"
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
    add_catalog_arguments(ap)
    ap.add_argument("--names", nargs="+", required=True, help="Member output ids")
    ap.add_argument("--members", nargs="+", required=True, help="Member artifact paths")
    ap.add_argument("--out-file", required=True, help="Bundle path (lc {output})")
    ap.add_argument("--harmonic-validation-binning", type=HARMONIC_VALIDATION_BINNING.__getitem__)
    a = ap.parse_args(argv)
    if a.layout == "cl" and a.harmonic_validation_binning is None:
        ap.error("--layout cl requires --harmonic-validation-binning")

    if len(a.names) != len(a.members):
        ap.error(
            f"--names ({len(a.names)}) and --members ({len(a.members)}) must be "
            "the same length; they come from the same expansion"
        )

    catalogs = catalogs_from(a)
    # Scan IDs describe fixed roles, even when the primary catalogue is raw.
    versions = catalogs.scan_slots
    expected = (_cl_map(versions, a.harmonic_validation_binning, catalogs)
                if a.layout == "cl" else LAYOUTS[a.layout](versions, catalogs))

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
