#!/usr/bin/env python3
"""Emit LC_TASK_CPUS_MAP from the sub-analysis specs' declared resources.

Login-node safe: stdlib only (pyyaml is not guaranteed present), so the
`resources: { cpus: N, ... }` lines are read with a regex keyed to the
enclosing `- id:`.  Only entries above one cpu are emitted, capped at the
allocation's 24 cpus per node.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

CPU_CAP = 24

_ID_RE = re.compile(r"^\s*-\s+id:\s*(\S+)\s*$")
_CPUS_RE = re.compile(r"^\s*resources:\s*\{[^}]*\bcpus:\s*(\d+)")


def cpus_by_output(spec: Path) -> dict[str, int]:
    """Declared cpus per output id, in the order the spec declares them."""

    found: dict[str, int] = {}
    current: str | None = None
    for line in spec.read_text().splitlines():
        id_match = _ID_RE.match(line)
        if id_match:
            current = id_match.group(1)
            continue
        cpus_match = _CPUS_RE.match(line)
        if cpus_match and current is not None:
            found[current] = int(cpus_match.group(1))
            current = None
    return found


def main(root: Path) -> int:
    mapping: dict[str, int] = {}
    for spec in sorted(root.glob("analyses/*/astra.yaml")):
        analysis = spec.parent.name
        for output_id, cpus in cpus_by_output(spec).items():
            if cpus > 1:
                mapping[f"{analysis}/{output_id}"] = min(cpus, CPU_CAP)
    json.dump(mapping, sys.stdout, sort_keys=True, separators=(",", ":"))
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()))
