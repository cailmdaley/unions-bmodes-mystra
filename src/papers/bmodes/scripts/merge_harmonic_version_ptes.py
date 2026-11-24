"""Gather the per-version harmonic COSEBI B-mode PTEs into one JSON.

`harmonic_config_cosebis_comparison.py` computes harmonic- and config-space
COSEBIs for one catalog version per output. This is the gather half of the `cosebis_harmonic_scan` astra-sweep: it takes the
per-version modes NPZs produced by the sweep members and emits, for both angular
ranges, the same `harmonic_b_mode_ptes` / `config_b_mode_ptes` dictionaries the
per-version metrics pass produces — merged across versions, so
`emit_cosebis_harmonic_config.py --harmonic-version-ptes` can fill the
non-fiducial harmonic rows.

    python merge_harmonic_version_ptes.py \
        --members <modes.npz>... --out-file <path.json>
"""

import argparse
import json
import os
import tempfile

from harmonic_config_cosebis_comparison import _load_modes, _metrics_from_modes

ANGULAR_RANGES = ("fiducial", "full")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--members", nargs="+", required=True, help="Per-version modes NPZ paths"
    )
    ap.add_argument("--out-file", required=True, help="Merged JSON path (lc {output})")
    a = ap.parse_args(argv)

    merged = {r: {"harmonic_b_mode_ptes": {}, "config_b_mode_ptes": {}} for r in ANGULAR_RANGES}
    versions = []

    for npz in a.members:
        modes = _load_modes(npz)
        version = str(modes["version"])
        versions.append(version)
        for angular_range in ANGULAR_RANGES:
            # _metrics_from_modes is the single source of truth for the PTE/chi2
            # convention (reliable-mode truncation, covariance block slicing);
            # go through it rather than re-deriving here.
            with tempfile.NamedTemporaryFile(
                "r+", suffix=".json", dir=os.environ.get("TMPDIR")
            ) as tmp:
                _metrics_from_modes(modes, angular_range, tmp.name)
                one = json.load(open(tmp.name))
            for key in ("harmonic_b_mode_ptes", "config_b_mode_ptes"):
                merged[angular_range][key].update(one[key])

    payload = dict(merged)
    payload["versions"] = versions
    with open(a.out_file, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Wrote {a.out_file} ({len(versions)} versions x {len(ANGULAR_RANGES)} cuts)")


if __name__ == "__main__":
    main()
