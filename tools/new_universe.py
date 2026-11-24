"""Add a universe that differs from an existing one in the decisions given.

A universe is one selection of every decision. This copies the parent's
(`baseline` unless --from says otherwise) at the root and in each sub-analysis,
changes the decisions named, and links Studio's results directories. Nothing is
computed: `lc materialize` does that.

    python tools/new_universe.py initial fiducial_catalog=initial
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main(name, parent, choices):
    specs = [ROOT, *(p.parent for p in sorted(ROOT.glob("analyses/**/astra.yaml")))]
    texts = {d / "universes" / f"{name}.yaml": (d / "universes" / f"{parent}.yaml").read_text()
             for d in specs if (d / "universes" / f"{parent}.yaml").exists()}
    for path, text in texts.items():
        text = re.sub(r"^id: .*$", f"id: {name}", text, flags=re.M)
        text = re.sub(r"^description: \|\n(?:  .*\n|\n)*", "", text, flags=re.M)
        chosen = ", ".join(f"{k} = {v}" for k, v in choices.items())
        text = text.replace(f"id: {name}\n", f"id: {name}\ndescription: |\n  {parent} with {chosen}.\n", 1)
        texts[path] = re.sub(rf"universe: {parent}$", f"universe: {name}", text, flags=re.M)
    for decision, option in choices.items():
        pattern = re.compile(rf"^  {decision}: .*$", re.M)
        hits = [path for path, text in texts.items() if pattern.search(text)]
        if not hits:
            sys.exit(f"no universe file selects {decision!r}")
        for path in hits:
            texts[path] = pattern.sub(f"  {decision}: {option}", texts[path])
    for path, text in texts.items():
        path.write_text(text)
    subprocess.run([sys.executable, ROOT / "tools/studio_links.py", name], check=True)
    print(f"universe {name}: {len(texts)} files")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("name")
    ap.add_argument("choices", nargs="+", metavar="decision=option")
    ap.add_argument("--from", dest="parent", default="baseline")
    a = ap.parse_args()
    main(a.name, a.parent, dict(c.split("=", 1) for c in a.choices))
