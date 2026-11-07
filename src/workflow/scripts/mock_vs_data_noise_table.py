"""Tabulate the mock-versus-data noise budget.

Each row reads the survey scale of one `reference_covariance` option from
`decision_values.REFERENCE_COVARIANCE`, which carries the area convention, the
galaxy count in density clothing (`n_eff`) and the total ellipticity dispersion,
and the nside-4096 footprint audit, which measures the two catalogues' occupied
sky and their raw row counts.

The derived columns are the whole argument of the density finding.  Shape
noise per galaxy is σ_e²/2 per component, so the noise *power* of a catalogue
scales as σ_e²/N, and a B-mode variance — which for these mocks is pure noise,
there being no B signal by construction — scales as the square of that.  The
area convention cancels from every ratio because both catalogues occupy one
sky, which the footprint audit is what establishes.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from decision_values import REFERENCE_COVARIANCE

ARCMIN2_PER_DEG2 = 3600.0


def read_footprint_audit(path):
    """Pull the two catalogues' row counts and occupied areas out of the audit."""

    text = Path(path).read_text()
    result = {}
    for label in ("data", "mock"):
        match = re.search(
            rf"^{label}:\s*N=(\d+),\s*occupied pix=\d+,\s*area=([0-9.]+) deg2",
            text,
            re.MULTILINE,
        )
        if match is None:
            raise ValueError(f"{path} has no '{label}:' galaxy/area line")
        result[f"{label}_galaxies"] = int(match.group(1))
        result[f"{label}_area_deg2"] = float(match.group(2))
    overlap = re.search(r"^overlap\s*=\s*([0-9.]+) deg2", text, re.MULTILINE)
    if overlap is None:
        raise ValueError(f"{path} has no 'overlap =' line")
    result["overlap_area_deg2"] = float(overlap.group(1))
    return result


def build_rows(data, mock, audit):
    """One row per catalogue, with the ratios taken against the data row."""

    def count(entry):
        # n_eff is a galaxy count wearing density units: the count divided by
        # the fixed area convention.  Multiplying it back out is the
        # convention-free reading.
        return entry.n_eff * entry.area * ARCMIN2_PER_DEG2

    def noise_power(entry):
        return (entry.sigma_e_combined / math.sqrt(2.0)) ** 2 / count(entry)

    reference = noise_power(data)
    rows = []
    for label, entry, galaxies, area in (
        (
            "UNIONS v1.4.6.3 leak-corrected catalogue",
            data,
            audit["data_galaxies"],
            audit["data_area_deg2"],
        ),
        (
            "GLASS v1.4.6.3_v2 mock suite",
            mock,
            audit["mock_galaxies"],
            audit["mock_area_deg2"],
        ),
    ):
        ratio = noise_power(entry) / reference
        rows.append(
            {
                "catalogue": label,
                "galaxies_counted": galaxies,
                "galaxies_from_n_gal_convention": round(count(entry)),
                "occupied_area_deg2": area,
                "area_convention_deg2": entry.area,
                "n_eff_arcmin2": entry.n_eff,
                "sigma_e_per_component": entry.sigma_e_combined / math.sqrt(2.0),
                "noise_power_ratio_to_data": ratio,
                "bb_variance_ratio_to_data": ratio**2,
            }
        )
    return rows


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--footprint-audit", required=True)
    parser.add_argument("--out-file", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = _parser(argv)
    audit = read_footprint_audit(args.footprint_audit)
    rows = build_rows(
        REFERENCE_COVARIANCE["data_covariance"],
        REFERENCE_COVARIANCE["mockmatched"],
        audit,
    )
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        f"{value:.6g}" if isinstance(value, float) else value
                    )
                    for key, value in row.items()
                }
            )
    print(
        f"Saved {out_path}: overlap {audit['overlap_area_deg2']} deg2, "
        f"noise-power ratio {rows[1]['noise_power_ratio_to_data']:.4f}, "
        f"BB variance ratio {rows[1]['bb_variance_ratio_to_data']:.4f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
