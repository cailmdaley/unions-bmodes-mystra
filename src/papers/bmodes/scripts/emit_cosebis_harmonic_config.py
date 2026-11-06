"""Emit a flat harmonic-vs-config COSEBI PTE comparison table.

Numbers quoted in the B-modes paper compare the COSEBIs B-mode PTE computed
two ways — harmonic space (pseudo-Cℓ bandpower fit) and config space
(real-space COSEBIs coefficients) — across catalog version and angular-range
cut. Those values live inside nested-dict `evidence.json` carriers belonging
to `type: figure` / `type: data` ASTRA outputs
(``null_tests.config_space_pte_evidence``, ``cosebis.cosebis_harmonic_modes``,
and ``cosebis.cosebis_harmonic_modes_full``). MySTRA's ``{astra:value col=
where=}`` role can only read a flat CSV/JSON table, not a nested dict — see
``work/transpile-map.md`` §0 in the ASTRA reproduction. This script reads
those *already-materialized* evidence.json files (it does not recompute
anything) and writes the requested cells out as a flat CSV table, one row per
(version, path, cut) triple.

Config-space PTEs are computed empirically (bootstrap resampling), so they
carry no chi2/dof — those columns are left blank for `path=config` rows.

    python emit_cosebis_harmonic_config.py \
        --config-space-evidence /path/to/config_space_pte_evidence/evidence.json \
        --cosebis-evidence /path/to/cosebis_harmonic_modes/evidence.json \
        --cosebis-evidence-full /path/to/cosebis_harmonic_modes_full/evidence.json \
        --out <output_dir>
"""

import argparse
import csv
import json
import math



def _sig3(x):
    """Format a float to 3 significant figures; scientific for |x|<1e-3 or >=1e4, else fixed."""
    if x == "" or x is None:
        return ""
    x = float(x)
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    if e < -3 or e >= 4:
        return f"{x:.2e}"  # mantissa carries 3 sig figs, e.g. 7.37e-10
    return f"{x:.{max(0, 2 - e)}f}"  # 3 sig figs fixed, e.g. 0.938, 0.00461


# Full leak-corrected version string -> CSV slug, in the paper's version order.
_VERSION_SLUGS = {
    "SP_v1.4.5_leak_corr": "initial",
    "SP_v1.4.6.3_leak_corr": "size_cuts",
    "SP_v1.4.8_leak_corr": "masked",
    "SP_v1.4.11.3_leak_corr": "relaxed_flags",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-space-evidence", required=True)
    parser.add_argument("--cosebis-evidence", required=True)
    parser.add_argument("--cosebis-evidence-full", required=True)
    parser.add_argument(
        "--harmonic-version-ptes",
        default=None,
        help=(
            "Merged per-version harmonic PTE JSON from the cosebis_harmonic_scan "
            "sweep gather ({'fiducial': {...}, 'full': {...}}); supplies the "
            "non-fiducial leak-corrected versions' harmonic rows."
        ),
    )
    parser.add_argument("--leakage-correction", required=True, choices=["applied", "uncorrected"])
    parser.add_argument("--out-file", required=True, help="output CSV path")
    args = parser.parse_args()

    def _payload(path):
        """Accept the legacy {'evidence': {...}} wrapper and the new
        top-level metrics format alike."""
        d = json.load(open(path))
        return d.get("evidence", d)

    config_evidence = _payload(args.config_space_evidence)
    harmonic_fid_evidence = _payload(args.cosebis_evidence)
    harmonic_full_evidence = _payload(args.cosebis_evidence_full)

    slugs = _VERSION_SLUGS
    if args.leakage_correction == "uncorrected":
        slugs = {k.removesuffix("_leak_corr"): v for k, v in slugs.items()}
    config_versions = config_evidence.get("versions", {})
    harmonic_fid_ptes = harmonic_fid_evidence.get("harmonic_b_mode_ptes", {})
    harmonic_full_ptes = harmonic_full_evidence.get("harmonic_b_mode_ptes", {})

    # The per-version metrics files carry the fiducial version only; the other
    # three leak-corrected versions come from the cosebis_harmonic_scan gather.
    # The fiducial evidence files win on collision.
    if args.harmonic_version_ptes:
        sweep = json.load(open(args.harmonic_version_ptes))
        harmonic_fid_ptes = {
            **sweep.get("fiducial", {}).get("harmonic_b_mode_ptes", {}),
            **harmonic_fid_ptes,
        }
        harmonic_full_ptes = {
            **sweep.get("full", {}).get("harmonic_b_mode_ptes", {}),
            **harmonic_full_ptes,
        }

    for version_full in slugs:
        if version_full not in harmonic_fid_ptes:
            print(
                f"WARNING: no harmonic PTEs for {version_full}; "
                "harmonic rows will be empty for this version "
                "(expected only if --harmonic-version-ptes was omitted)."
            )

    rows = []
    for version_full, slug in slugs.items():
        # Harmonic space: fiducial and full-range come from separate evidence
        # carriers (--angular-range fiducial vs full).
        for cut, harm_ptes in [
            ("fiducial", harmonic_fid_ptes),
            ("full_range", harmonic_full_ptes),
        ]:
            harm = harm_ptes.get(version_full, {})
            rows.append(
                {
                    "version": slug,
                    "path": "harmonic",
                    "cut": cut,
                    "pte": _sig3(harm.get("pte")),
                    "chi2": _sig3(harm.get("chi2")),
                    "dof": harm.get("dof"),
                }
            )

        # Config space: single evidence carrier, both cuts as separate keys.
        cfg = config_versions.get(version_full, {}).get("cosebis_stats", {})
        for cut, pte_key in [
            ("fiducial", "pte_at_fiducial"),
            ("full_range", "pte_at_full_range"),
        ]:
            rows.append(
                {
                    "version": slug,
                    "path": "config",
                    "cut": cut,
                    "pte": _sig3(cfg.get(pte_key)),
                    "chi2": "",
                    "dof": "",
                }
            )

    fieldnames = ["version", "path", "cut", "pte", "chi2", "dof"]
    with open(args.out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.out_file} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
