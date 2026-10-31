"""Tabulate the mock-versus-data noise budget from the artifacts that fix it.

Every number here is read or derived from three files and nothing else:

* the two CosmoCov ini files — the data catalogue's survey entry and the
  `SP_v1.4.6.3_mockv2` entry generated at the mocks' own noise — which carry
  the area convention, the galaxy count in density clothing (`source_n_gal`)
  and the total ellipticity dispersion (`sigma_e`);
* the nside-4096 footprint audit, which measures the two catalogues' occupied
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

ARCMIN2_PER_DEG2 = 3600.0

_INI_KEYS = ("area", "sigma_e", "source_n_gal")


def read_cosmocov_ini(path):
    """Pull the survey block's area, sigma_e and source_n_gal out of an ini."""

    values = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"\s*([A-Za-z_]+)\s*:\s*([0-9.eE+-]+)\s*$", line)
        if match and match.group(1) in _INI_KEYS and match.group(1) not in values:
            values[match.group(1)] = float(match.group(2))
    missing = [key for key in _INI_KEYS if key not in values]
    if missing:
        raise KeyError(f"{path} declares no {missing}")
    return values


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


def build_rows(data_ini, mock_ini, audit):
    """One row per catalogue, with the ratios taken against the data row."""

    def count(entry):
        # source_n_gal is a galaxy count wearing density units: the count
        # divided by the fixed area convention.  Multiplying it back out is
        # the convention-free reading.
        return entry["source_n_gal"] * entry["area"] * ARCMIN2_PER_DEG2

    def noise_power(entry):
        return (entry["sigma_e"] / math.sqrt(2.0)) ** 2 / count(entry)

    reference = noise_power(data_ini)
    rows = []
    for label, entry, galaxies, area in (
        (
            "UNIONS v1.4.6.3 leak-corrected catalogue",
            data_ini,
            audit["data_galaxies"],
            audit["data_area_deg2"],
        ),
        (
            "GLASS v1.4.6.3_v2 mock suite",
            mock_ini,
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
                "area_convention_deg2": entry["area"],
                "n_eff_arcmin2": entry["source_n_gal"],
                "sigma_e_per_component": entry["sigma_e"] / math.sqrt(2.0),
                "noise_power_ratio_to_data": ratio,
                "bb_variance_ratio_to_data": ratio**2,
            }
        )
    return rows


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cosmocov-ini-data", required=True)
    parser.add_argument("--cosmocov-ini-mock", required=True)
    parser.add_argument("--footprint-audit", required=True)
    parser.add_argument("--out-file", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = _parser(argv)
    audit = read_footprint_audit(args.footprint_audit)
    rows = build_rows(
        read_cosmocov_ini(args.cosmocov_ini_data),
        read_cosmocov_ini(args.cosmocov_ini_mock),
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
