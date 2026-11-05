"""Check local Markdown links in the report against its explicit MyST labels.

The publication theme provides the ASTRA inventory target. Scientific equations
also declare labels through LaTeX or MyST's trailing equation-label syntax.

    uv run python tools/check_report_links.py [index.md]
"""
import re
import sys
from pathlib import Path


def broken_links(text):
    labels = set(re.findall(r"^\(([^\s)]+)(?:\s+[^)]*)?\)=", text, re.M))
    labels.update(re.findall(r"\\label\{([^}]+)\}", text))
    labels.update(re.findall(r"^\$\$\s+\(([^)]+)\)", text, re.M))
    labels.add("astra-inventory")
    return sorted(set(re.findall(r"\]\(#([^\s)]+)\)", text)) - labels)


if __name__ == "__main__":
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "index.md")
    missing = broken_links(path.read_text())
    for target in missing:
        print(f"{path}: unresolved report link #{target}")
    sys.exit(bool(missing))
