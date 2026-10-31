"""The fiducial catalog's products, resolved from its Paper II config.

The catalog config (`config/catalogs/<fiducial_catalog>.yaml`) names the
cat_config version (`fiducial.version`, leak-corrected) and the footprint mask
power spectrum (`fiducial.mask_cls`); cat_config supplies the rest, with the
n(z) file switched to `--blind` the way sp_validation does (`_<blind>.txt`).
One product per call, written to `--out-file`:

    redshift_distribution  copy of the n(z) text file
    mask_power_spectrum    copy of the footprint mask power spectrum
    survey_properties      JSON {area, n_eff, sigma_e_combined} from cov_th
"""

import argparse
import json
import re
import shutil

import yaml

from paper_config import load_config

EMITS = ("redshift_distribution", "mask_power_spectrum", "survey_properties")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, help="Catalog config (config/catalogs/<option>.yaml)")
    ap.add_argument("--cat-config", required=True)
    ap.add_argument("--emit", required=True, choices=EMITS)
    ap.add_argument(
        "--blind",
        type=lambda s: s.removeprefix("blind_").upper(),
        choices=["A", "B", "C"],
        default="A",
        help="A/B/C, or the fiducial_blind option (blind_a/blind_b/blind_c)",
    )
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    fiducial = load_config(a.config)["fiducial"]
    with open(a.cat_config) as f:
        entry = yaml.safe_load(f)[fiducial["version"].removesuffix("_leak_corr")]

    if a.emit == "redshift_distribution":
        shutil.copy(re.sub(r"_[ABC]\.txt$", f"_{a.blind}.txt", entry["shear"]["redshift_path"]), a.out_file)
    elif a.emit == "mask_power_spectrum":
        shutil.copy(fiducial["mask_cls"], a.out_file)
    else:
        c = entry["cov_th"]
        with open(a.out_file, "w") as f:
            json.dump(dict(area=c["A"], n_eff=c["n_e"], sigma_e_combined=c["sigma_e"]), f, indent=2)


if __name__ == "__main__":
    main()
