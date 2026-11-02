"""Check that a spec diff touches only prose fields.

Compares every hand-written astra.yaml against a git revision with the prose
fields (description, rationale, label, claim, notes, scope, excluded_reason)
removed. Anything else that differs — ids, recipes, inputs, decisions,
evidence quotes, defaults — is reported and the exit status is 1.

Findings may be renamed or replaced when --allow-findings is passed.

    python tools/check_prose_only.py [REV] [--allow-findings]
"""
import subprocess
import sys
from pathlib import Path

import yaml

PROSE = {"description", "rationale", "label", "claim", "notes", "scope", "excluded_reason"}


def strip(node, allow_findings):
    if isinstance(node, dict):
        return {
            k: strip(v, allow_findings)
            for k, v in node.items()
            if k not in PROSE and not (allow_findings and k == "findings")
        }
    if isinstance(node, list):
        return [strip(v, allow_findings) for v in node]
    return node


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    rev = args[0] if args else "HEAD"
    allow = "--allow-findings" in sys.argv
    files = subprocess.run(
        ["git", "ls-files", "astra.yaml", "analyses/*/astra.yaml"],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    bad = 0
    for f in files:
        old = subprocess.run(["git", "show", f"{rev}:./{f}"], capture_output=True, text=True)
        if old.returncode:
            print(f"new file, not compared: {f}")
            continue
        a = strip(yaml.safe_load(old.stdout), allow)
        b = strip(yaml.safe_load(Path(f).read_text()), allow)
        if a != b:
            bad += 1
            print(f"non-prose change: {f}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
