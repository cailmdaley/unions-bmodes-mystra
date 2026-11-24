#!/usr/bin/env python3
"""Convert generated paper macros to CSV, or compare two such tables."""

import argparse
import csv
import re
from pathlib import Path


NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
MACRO = re.compile(r"\\newcommand\{\\([A-Za-z]+)\}\{")
SUPERSCRIPT = str.maketrans("-+0123456789", "⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹")


def definitions(text):
    for match in MACRO.finditer(text):
        depth = 1
        start = match.end()
        i = start
        while i < len(text) and depth:
            if text[i] == "{" and (i == 0 or text[i - 1] != "\\"):
                depth += 1
            elif text[i] == "}" and (i == 0 or text[i - 1] != "\\"):
                depth -= 1
            i += 1
        if depth:
            raise ValueError(f"Unclosed macro {match.group(1)}")
        yield match.group(1), text[start:i - 1]


def render(body):
    number = re.fullmatch(r"\\num\{(" + NUMBER + r")\}", body)
    if number:
        raw = number.group(1)
        parts = re.split("[eE]", raw)
        display = raw if len(parts) == 1 else f"{parts[0]} × 10{str(int(parts[1])).translate(SUPERSCRIPT)}"
        return display, float(raw)
    percent = re.fullmatch(r"(" + NUMBER + r")\\%", body)
    if percent:
        return percent.group(1) + "%", float(percent.group(1))
    if re.fullmatch(NUMBER, body):
        return body, float(body)
    # Preserve math as MyST; unlike numerical claims it has no scalar value.
    return body.replace(r"\%", "%"), ""


def write_rows(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def extract(source, destination):
    rows = []
    seen = set()
    for name, body in definitions(Path(source).read_text(encoding="utf-8")):
        if name in seen:
            raise ValueError(f"Duplicate macro: {name}")
        seen.add(name)
        display, value = render(body)
        rows.append(dict(macro=name, display=display, value=value))
    if not rows:
        raise ValueError("No macro definitions found")
    write_rows(destination, ["macro", "display", "value"], rows)


def compare(current, submitted, destination):
    def read(path):
        with open(path, newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        result = {row["macro"]: row["display"] for row in rows}
        if len(result) != len(rows):
            raise ValueError(f"Duplicate macro in {path}")
        return result
    accepted, reference = read(current), read(submitted)
    rows = [dict(macro=name, submitted=reference.get(name, ""),
                 accepted=accepted.get(name, ""),
                 changed=str(reference.get(name) != accepted.get(name)).lower())
            for name in sorted(accepted.keys() | reference.keys())]
    write_rows(destination, ["macro", "submitted", "accepted", "changed"], rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    extract_parser = commands.add_parser("extract")
    extract_parser.add_argument("source")
    extract_parser.add_argument("destination")
    compare_parser = commands.add_parser("compare")
    compare_parser.add_argument("current")
    compare_parser.add_argument("submitted")
    compare_parser.add_argument("destination")
    args = parser.parse_args()
    if args.command == "extract":
        extract(args.source, args.destination)
    else:
        compare(args.current, args.submitted, args.destination)


if __name__ == "__main__":
    main()
