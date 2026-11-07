"""Check that every spec cross-link `[text](#...kind.id)` resolves.

A link is read from the spec holding it: each `../` climbs one analysis,
`analyses.<id>.` or `<id>.` descends into a sub-analysis, and the tail names an
element that analysis declares (`kind.id`). A link to a whole sub-analysis is
refused, because report dialogs have nothing to open for it. Prose that a scan's
sidecar generates is checked in the generated astra.yaml.

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


def resolves(spec_dir, target):
    """Whether `target` names something, read from the spec in `spec_dir`."""
    here = spec_dir
    while target.startswith("../"):
        here, target = here.parent.parent if here != Path(".") else here, target[3:]
    parts = target.split(".")
    while parts:
        if parts[0] in KINDS:
            if parts[0] == "analyses" or len(parts) != 2:
                return False
            return parts[1] in ids(yaml.safe_load((here / "astra.yaml").read_text()), parts[0])
        child = parts[1] if parts[0] == "analyses" and len(parts) > 1 else parts[0]
        parts = parts[2:] if parts[0] == "analyses" else parts[1:]
        here = here / "analyses" / child
        if not (here / "astra.yaml").exists():
            return False
    return False


def main():
    bad = 0
    for spec in [Path("astra.yaml"), *sorted(Path("analyses").rglob("astra.yaml"))]:
        text = re.sub(r"\s*\n\s*", " ", "\n".join(
            l for l in spec.read_text().splitlines() if not l.lstrip().startswith("#")
        ))
        for link in LINK.finditer(text):
            if not resolves(spec.parent, link.group(1)):
                print(f"{spec}: dead link #{link.group(1)}")
                bad += 1
    print(f"links: {bad} dead")
    sys.exit(bool(bad))


if __name__ == "__main__":
    main()
