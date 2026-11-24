"""Verify that MyST preserves the report's bound values in its built AST.

Run after `myst build`:
    uv run --no-project python tools/check_report_build.py

ASTRA roles inside caption overrides aren't resolved by MySTRA v0.0.8.
This check detects lost bindings even when the Markdown source is unchanged.
"""
import json
import re
import sys
from pathlib import Path


def walk(node):
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from walk(value)
    elif isinstance(node, list):
        for value in node:
            yield from walk(value)


def check(text, document):
    expected = len(re.findall(r"\{astra:value\b[^}]*\}`[^`]+`", text))
    nodes = list(walk(document))
    resolved = sum(
        node.get("data", {}).get("astra", {}).get("kind") == "value"
        for node in nodes
        if isinstance(node.get("data"), dict)
    )
    unresolved = sum(
        node.get("type") == "mystRole" and node.get("name", "").startswith("astra")
        for node in nodes
    )
    errors = []
    if resolved != expected:
        errors.append(f"bound values: {resolved} resolved, {expected} authored")
    if unresolved:
        errors.append(f"unresolved ASTRA roles: {unresolved}")
    return errors


if __name__ == "__main__":
    errors = check(
        Path("index.md").read_text(),
        json.loads(Path("_build/site/content/index.json").read_text()),
    )
    for error in errors:
        print(error)
    sys.exit(bool(errors))
