"""Check that every spec cross-link `[text](#...kind.id)` names an element that exists.

kind is one of inputs, outputs, decisions, findings, prior_insights; the id must
be declared under that kind in some hand-written or generated astra.yaml.

    python tools/check_links.py
"""
import re
import sys
from pathlib import Path

import yaml

KINDS = ("inputs", "outputs", "decisions", "findings", "prior_insights")
LINK = re.compile(r"\]\(\s*#([^)\s]+)\)")


def ids(spec, kind):
    block = spec.get(kind) or []
    if isinstance(block, dict):
        return set(block)
    return {e["id"] for e in block if isinstance(e, dict) and "id" in e}


def main():
    files = [Path("astra.yaml"), *sorted(Path("analyses").rglob("astra.yaml"))]
    known = {k: set() for k in KINDS}
    for f in files:
        spec = yaml.safe_load(f.read_text())
        for k in KINDS:
            known[k] |= ids(spec, k)
    bad = 0
    for f in files:
        text = f.read_text()
        for m in LINK.finditer(re.sub(r"\s*\n\s*", " ", text)):
            parts = m.group(1).split(".")
            for i, p in enumerate(parts[:-1]):
                if p.lstrip("./") in KINDS and parts[i + 1] not in known[p.lstrip("./")]:
                    print(f"{f}: dead link #{m.group(1)}")
                    bad += 1
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
