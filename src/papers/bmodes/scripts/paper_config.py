"""Load a Paper II config, resolving its `base:` chain.

Each `fiducial_catalog` option names a config under `config/catalogs/`
(`catalogs/<option>.yaml`) that inherits the shared `config.yaml` through
`base:` and overrides what the catalog changes: the version set, the fiducial
version, and the labels. Mappings merge key by key; any other value, lists
included, replaces the inherited one, and so does a mapping tagged `!replace`.

As a CLI it prints one dotted key, which is how recipes resolve the catalog's
cat_config version:

    python paper_config.py --config config/catalogs/size_cut.yaml --get fiducial.version
"""

import argparse
from pathlib import Path

import yaml


class _Replace(dict):
    """A mapping that replaces the inherited one instead of merging into it."""


class _Loader(yaml.SafeLoader):
    pass


_Loader.add_constructor(
    "!replace", lambda loader, node: _Replace(loader.construct_mapping(node, deep=True))
)


def _merge(base, override):
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, _Replace):
            merged[key] = dict(value)
        elif isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path):
    """The config at `path`, merged over its `base:` chain (paths relative to it)."""
    path = Path(path)
    with open(path) as f:
        config = yaml.load(f, Loader=_Loader) or {}
    base = config.pop("base", None)
    if base is None:
        return config
    return _merge(load_config(path.parent / base), config)


def select_leakage(config, correction):
    """Select the primary shear catalogue, retaining both variants in scans."""
    if correction == "applied":
        return config
    config["fiducial"]["version"] = config["fiducial"]["version"].removesuffix("_leak_corr")
    config["plotting"]["fiducial_for_comparison"] = config["fiducial"]["version"]
    config["plotting"]["version_labels"] = {
        k.removesuffix("_leak_corr"): v
        for k, v in config["plotting"]["version_labels"].items()
    }
    config["leakage_correction"] = correction
    return config


def scan_version(config, version, fine=False):
    """Resolve a stable scan slot; the primary variant is produced separately."""
    if config.get("leakage_correction") != "uncorrected":
        return version
    if fine:
        return version.removesuffix("_leak_corr")
    fiducial = config["fiducial"]["version"]
    return fiducial + "_leak_corr" if version == fiducial else version


def comparison_versions(config):
    corrected = config.get("leakage_correction") != "uncorrected"
    return [v for v in config["versions"]
            if "_ecut" not in v and ("_leak_corr" in v) == corrected]


def get(config, dotted):
    value = config
    for key in dotted.split("."):
        value = value[key]
    return value


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(description="Print one key of a Paper II config.")
    ap.add_argument("--config", required=True)
    ap.add_argument("--get", help="Dotted key, e.g. fiducial.version")
    ap.add_argument("--leakage-correction", choices=["applied", "uncorrected"])
    ap.add_argument("--out-file", help="Write the resolved analysis config")
    ap.add_argument("--scan-version", help="Resolve a fixed scan slot's catalogue")
    ap.add_argument("--fine-grid", action="store_true")
    ap.add_argument(
        "--base-version",
        action="store_true",
        help="Strip `_leak_corr`: the cat_config entry a leak-corrected version reads",
    )
    a = ap.parse_args(argv)
    config = load_config(a.config)
    if a.leakage_correction:
        config = select_leakage(config, a.leakage_correction)
    if a.out_file:
        Path(a.out_file).write_text(yaml.safe_dump(config, sort_keys=False))
        return
    value = scan_version(config, a.scan_version, a.fine_grid) if a.scan_version else get(config, a.get)
    if a.base_version:
        value = str(value).removesuffix("_leak_corr")
    print(value)


if __name__ == "__main__":
    _from_cli()
