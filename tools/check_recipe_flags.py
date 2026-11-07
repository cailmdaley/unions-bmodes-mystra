"""Check that every recipe calls its scripts with flags they accept.

Each `python <script>.py` and `bash <script>.sh` call in a recipe is compared
with the script's own parser: `add_argument` in Python, `--flag)` case arms in
shell. A flag the script does not define fails, and so does a required flag the
recipe leaves out. Command substitutions are checked as calls of their own.

    python tools/check_recipe_flags.py
"""
import ast
import re
import sys
from pathlib import Path

import yaml

CALL = re.compile(r"\b(?:python|bash)\s+(src/\S+\.(?:py|sh))")
FLAG = re.compile(r"(?<![\w-])--[a-z0-9][a-z0-9-]*")
CASE_ARM = re.compile(r"^\s*((?:--[a-z0-9-]+\|?)+)\)", re.M)
SCRIPT_DIRS = [Path("src/papers/bmodes/scripts"), Path("src/workflow/scripts")]


def parser_arguments(tree):
    """Yield (flag names, keywords, in_helper) for every `add_argument` in a module."""
    helpers = {
        id(call)
        for function in ast.walk(tree) if isinstance(function, ast.FunctionDef)
        if function.args.args and "parser" in function.args.args[0].arg
        for call in ast.walk(function)
    }
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "add_argument":
            names = [a.value for a in node.args if isinstance(a, ast.Constant) and str(a.value).startswith("--")]
            yield names, {k.arg: k.value for k in node.keywords}, id(node) in helpers


def local_imports(tree):
    """Sibling script modules a module imports, as paths."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return [d / f"{n}.py" for n in names for d in SCRIPT_DIRS if (d / f"{n}.py").exists()]


def python_flags(path):
    """Flags a Python script accepts, its own and those added by shared parser helpers."""
    tree = ast.parse(path.read_text())
    flags, required = set(), set()
    for names, keywords, in_helper in parser_arguments(tree):
        flags.update(names)
        if names and not in_helper and getattr(keywords.get("required"), "value", False):
            required.add(names[0])
    for module in local_imports(tree):
        for names, _, in_helper in parser_arguments(ast.parse(module.read_text())):
            if in_helper:
                flags.update(names)
    return flags, required


def shell_flags(path):
    """Flags a shell script parses, or None when it forwards `$@` to another program."""
    text = path.read_text()
    flags = {f for arm in CASE_ARM.findall(text) for f in arm.split("|")}
    return flags or None, set()


def calls(command):
    """Each script call with the text of its arguments, innermost substitutions first."""
    pieces = []
    while (found := re.findall(r"\$\(([^()]*)\)", command)):
        pieces += found
        command = re.sub(r"\$\(([^()]*)\)", "", command)
    for text in [*pieces, command]:
        for segment in re.split(r"\s(?:&&|\|\||;|\|)\s", text):
            match = CALL.search(segment)
            if match:
                yield match[1], segment[match.end():]


def main():
    bad = 0
    seen = set()
    for spec in [Path("astra.yaml"), *sorted(Path("analyses").rglob("astra.yaml"))]:
        for output in yaml.safe_load(spec.read_text()).get("outputs") or []:
            command = (output.get("recipe") or {}).get("command")
            if not command:
                continue
            command = "\n".join(l for l in command.splitlines() if not l.lstrip().startswith("#"))
            for line in command.replace("\\\n", " ").splitlines():
                for script, args in calls(" ".join(line.split())):
                    path = Path(script)
                    if not path.exists():
                        print(f"{spec}: {output['id']}: no such script {script}")
                        bad += 1
                        continue
                    flags, required = python_flags(path) if script.endswith(".py") else shell_flags(path)
                    if flags is None or (not flags and script.endswith(".py")):
                        continue
                    used = set(FLAG.findall(args))
                    problems = sorted(used - flags) + [f"missing {f}" for f in sorted(required - used)]
                    key = (script, tuple(problems))
                    if problems and key not in seen:
                        seen.add(key)
                        print(f"{spec}: {output['id']}: {script}: {', '.join(problems)}")
                        bad += 1
    print(f"recipe flags: {bad} problem(s)")
    sys.exit(bool(bad))


if __name__ == "__main__":
    main()
