"""Emit the COSEBI first-mode B-mode significance (B_1 / sigma_1) as an
ASTRA `type: metric` output.

The paper reports that the fiducial catalog's COSEBI B-mode first mode (n=1)
shows a ">4 sigma" excess on the full angular range (1-250 arcmin), dropping
below 1 sigma at the fiducial scale cuts. That significance is B_1/sigma_1,
where B_1 is the first COSEBI B-mode amplitude and sigma_1 its uncertainty.
Both live in the config-space COSEBI modes npz produced by
`cosebis_data_vector.py` (the source of the `fig_cosebis_fiducial` figure).

This is a thin transform: it reads the already-materialized npz and writes a
metric JSON (+ provenance evidence.json). No recomputation of the COSEBIs.
"""

import argparse
import json
import os

import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cosebis-modes-npz",
        required=True,
        help="cosebis_modes_<version>.npz (config-space COSEBIs)",
    )
    p.add_argument("--out-file", required=True, help="output metric JSON path")
    args = p.parse_args()

    d = np.load(args.cosebis_modes_npz)
    full_sig = float(d["full_Bn"][0] / d["full_sigma_B"][0])
    fid_sig = float(d["fiducial_Bn"][0] / d["fiducial_sigma_B"][0])
    version = str(d["version"]) if "version" in d else "unknown"

    # Metric artifact MySTRA reads: {value, uncertainty, unit}. `value` is the
    # headline full-range first-mode significance quoted in the figure caption;
    # provenance fields ride alongside (extra keys are ignored by readers).
    metric = {
        "value": round(full_sig, 3),
        "uncertainty": None,
        "unit": "sigma",
        "version": version,
        "full_range_significance": full_sig,
        "fiducial_cut_significance": fid_sig,
        "mode": 1,
    }
    with open(args.out_file, "w") as f:
        json.dump(metric, f, indent=2)

    print(f"full_range B1/sigma1 = {full_sig:.4f}  fiducial = {fid_sig:.4f}")


if __name__ == "__main__":
    main()
