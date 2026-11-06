"""The fiducial catalog's products, read from its cat_config entry.

`fiducial_catalog` names the cat_config version. Its entry gives the footprint
mask power spectrum (`mask_cls`) and the rest, with the n(z) file switched to
`--blind` the way sp_validation does (`_<blind>.txt`).
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

from decision_values import FIDUCIAL_BLIND, FIDUCIAL_CATALOG

EMITS = ("redshift_distribution", "mask_power_spectrum", "survey_properties")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fiducial-catalog", type=FIDUCIAL_CATALOG.__getitem__, required=True)
    ap.add_argument("--cat-config", required=True)
    ap.add_argument("--emit", required=True, choices=EMITS)
    ap.add_argument("--blind", type=FIDUCIAL_BLIND.__getitem__)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    with open(a.cat_config) as f:
        entry = yaml.safe_load(f)[a.fiducial_catalog]

    if a.emit == "redshift_distribution":
        shutil.copy(re.sub(r"_[ABC]\.txt$", f"_{a.blind}.txt", entry["shear"]["redshift_path"]), a.out_file)
    elif a.emit == "mask_power_spectrum":
        shutil.copy(entry["mask_cls"], a.out_file)
    else:
        c = entry["cov_th"]
        with open(a.out_file, "w") as f:
            json.dump(dict(area=c["A"], n_eff=c["n_e"], sigma_e_combined=c["sigma_e"]), f, indent=2)


if __name__ == "__main__":
    main()
