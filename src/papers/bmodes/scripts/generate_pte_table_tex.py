"""Render the Paper II PTE table as a LaTeX tabular."""

import argparse
import csv
import math
from pathlib import Path

_S_COLUMN = r"S[table-format=1.1e-2]"


def _format_two_sig(value):
    """Format a PTE with two significant figures, as input to a siunitx S column.

    PTEs at or above 0.1 are written in fixed notation (0.xx), smaller ones in
    scientific notation (a.b x 10^-n). The S column reserves one decimal: a
    fixed value's second decimal sits in the exponent space it leaves empty,
    so mantissas and exponents both align on the decimal marker.
    """
    value = float(value)
    if value != value:
        return "{--}"
    if value == 0:
        return "0"
    if abs(value) < 0.1:
        return f"{value:.1e}"
    decimals = max(0, 1 - int(math.floor(math.log10(abs(value)))))
    return f"{value:.{decimals}f}"


_VERSION_LABELS = {
    "Initial": "Initial",
    "Size_cuts": "Size cuts (fiducial)",
    "Masked": "Masked",
    "Relaxed_flags": "Relaxed flags",
}
_STAT_COLUMNS = [
    "cosebis_n6",
    "cosebis_n20",
    "xip_B",
    "xim_B",
    "xitot_B",
    "cl_BB",
]


def render_table(csv_path: Path, out_file: Path) -> None:
    with csv_path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    lines = [
        "% Auto-generated mock-calibrated PTE comparison tabular",
        "% Regenerate: generate_pte_table_tex.py",
        "% Wrap in table*/table environment in main tex for float control",
        r"\begin{tabular}{ll "
        + 2 * _S_COLUMN
        + r" @{\hskip 8pt} "
        + 3 * _S_COLUMN
        + r" @{\hskip 8pt} "
        + 2 * _S_COLUMN
        + "}",
        r"    \hline\hline",
        r"    \noalign{\vskip 3pt}",
        "    & & \\multicolumn{2}{c}{COSEBIs} & \\multicolumn{3}{c}{Pure E/B} & {$C_\\ell$} & \\\\",
        r"    \cmidrule(lr){3-4} \cmidrule(lr){5-7} \cmidrule(l){8-8} \cmidrule(l){9-9}",
        "    Version & Scale cuts & {$B_n$ ($n \\leq 6$)} & {$B_n$ ($n \\leq 20$)} & {$\\xi_+^{\\mathrm{B}}$} & {$\\xi_-^{\\mathrm{B}}$} & {$\\xi_{\\mathrm{tot}}^{\\mathrm{B}}$} & {$C_\\ell^{BB}$} & {min PTE} \\\\",
        r"    \hline",
        r"    \noalign{\vskip 3pt}",
    ]

    previous_version = None
    for row in rows:
        version = row["version"]
        version_label = _VERSION_LABELS.get(version, version.replace("_", " "))
        version_cell = version_label if version != previous_version else ""
        previous_version = version
        cut_label = row["cut"].replace("_", " ")
        values = [float(row[column]) for column in _STAT_COLUMNS]
        formatted = [_format_two_sig(value) for value in values]
        min_pte = min(values)
        formatted.append(_format_two_sig(min_pte))
        lines.append(
            "    " + " & ".join([version_cell, cut_label, *formatted]) + r"  \\"
        )
        if cut_label == "Full range":
            lines.append(r"    \noalign{\vskip 2pt}")

    if lines[-1].startswith(r"    \noalign"):
        lines.pop()
    lines.extend([r"  \hline", r"\end{tabular}"])
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"  → {out_file}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pte-csv", required=True)
    parser.add_argument("--out-file", required=True)
    args = parser.parse_args(argv)
    render_table(Path(args.pte_csv), Path(args.out_file))


if __name__ == "__main__":
    main()
