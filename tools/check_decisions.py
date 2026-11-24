"""Check that every decision reaches a reported result and tags form families.

A reported result is a root output whose type, following `from` aliases, is not
data. A decision reaches one when the result's dependency closure (`inputs`, an
input's `from`, an output's re-export `from`) contains an output declaring it.
A decision with `from` is the decision it aliases and carries no tags. Every
other decision carries exactly one family tag, shared with at least one other
decision and distinct from every analysis id.

    python tools/check_decisions.py
"""
import sys
from collections import Counter
from pathlib import Path

from check_reachable import deps, elems, load, specs, up


def canonical(scope, id_):
    elem = elems(scope, "decisions")[id_]
    return canonical(*up(scope, elem["from"])) if "from" in elem else (scope, id_)


def source(node):
    elem = elems(node[0], node[1])[node[2]]
    return source(deps(node)[0]) if "from" in elem else elem


def reached():
    stack = [
        ((), "outputs", oid)
        for oid in elems((), "outputs")
        if source(((), "outputs", oid)).get("type") != "data"
    ]
    seen, found = set(), set()
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack += deps(node)
        if node[1] == "outputs":
            for d in elems(node[0], "outputs")[node[2]].get("decisions") or []:
                found.add(canonical(node[0], d))
    return found


def problems():
    decisions = {
        (scope, id_): elem
        for scope in specs
        for id_, elem in elems(scope, "decisions").items()
        if "from" not in elem
    }
    name = lambda key: ".".join(key[0] + (key[1],))
    hit = reached()
    for key in decisions:
        if key not in hit:
            yield f"decision reaches no reported result: {name(key)}"
    for key, elem in decisions.items():
        if len(elem.get("tags") or []) != 1:
            yield f"decision needs exactly one family tag: {name(key)}"
    counts = Counter(t for elem in decisions.values() for t in elem.get("tags") or [])
    analyses = {scope[-1] for scope in specs if scope}
    for tag, n in sorted(counts.items()):
        if n == 1:
            yield f"tag used by one decision: {tag}"
        if tag in analyses:
            yield f"tag equals an analysis id: {tag}"


def main():
    load((), Path("astra.yaml"))
    found = list(problems())
    print("\n".join(found) if found else f"decisions ok: {len(specs)} specs")
    sys.exit(bool(found))


if __name__ == "__main__":
    main()
