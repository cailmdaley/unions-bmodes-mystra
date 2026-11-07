"""Check that every declared output is used by the paper.

An output is used when the report (index.md) renders or cites it, a finding
cites it as evidence, it is listed in KEEP below, or a used output consumes it
(through `inputs`, an input's `from`, or an output's re-export `from`).
Prints every unused output and exits 1 if there is one.

    python tools/check_reachable.py
"""
import re
import sys
from pathlib import Path

import yaml

# Reported outside the report's rendered elements.
KEEP = {
    (): {
        "pte_table_tex",  # the paper's TeX appendix table
        "primary_null_test_report",  # the primary-report selection the README documents
        "bb_covariance_blind_independence_evidence",  # §3.7 blind-independence statement
        "tenk_grid_comparison",  # evidence for the integration-grid convergence findings
    },
}

specs = {}


def load(scope, path):
    spec = yaml.safe_load(path.read_text())
    specs[scope] = spec
    for sub in spec.get("analyses") or {}:
        child = path.parent / "analyses" / sub / "astra.yaml"
        if child.exists():
            load(scope + (sub,), child)


def elems(scope, kind):
    block = specs[scope].get(kind) or []
    return {e["id"]: e for e in block} if isinstance(block, list) else block


def up(scope, ref):
    while ref.startswith("../"):
        scope, ref = scope[:-1], ref[3:]
    return scope, ref


def resolve(scope, ref):
    if "." in ref:
        parts = ref.split(".")
        return [(scope + tuple(parts[:-1]), "outputs", parts[-1])]
    for kind in ("outputs", "inputs"):
        if ref in elems(scope, kind):
            return [(scope, kind, ref)]
    return []


def deps(node):
    scope, kind, id_ = node
    elem = elems(scope, kind).get(id_) if scope in specs else None
    if elem is None:
        return []
    out = []
    if "from" in elem:
        src = elem["from"]
        out += resolve(*(up(scope, src) if src.startswith("../") else (scope, src)))
    for ref in elem.get("inputs") or []:
        out += resolve(scope, ref)
    return out


def main():
    load((), Path("astra.yaml"))
    roots = set()
    for m in re.finditer(r"(?:analyses\.([a-z_]+)\.)?outputs\.([a-z0-9_]+)", Path("index.md").read_text()):
        roots.add(((m.group(1),) if m.group(1) else (), "outputs", m.group(2)))
    for scope in specs:
        for finding in elems(scope, "findings").values():
            for ev in finding.get("evidence") or []:
                if "artifact" in ev:
                    roots |= set(resolve(scope, ev["artifact"]))
    for scope, ids in KEEP.items():
        roots |= {(scope, "outputs", i) for i in ids}
    seen, stack = set(), list(roots)
    while stack:
        node = stack.pop()
        if node not in seen:
            seen.add(node)
            stack += deps(node)
    unused = [
        ".".join(scope + (oid,))
        for scope in specs
        for oid in elems(scope, "outputs")
        if (scope, "outputs", oid) not in seen
    ]
    for u in unused:
        print(f"unused output: {u}")
    sys.exit(1 if unused else 0)


if __name__ == "__main__":
    main()
