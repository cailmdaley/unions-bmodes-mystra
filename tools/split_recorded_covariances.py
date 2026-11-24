"""Split a recorded pure-mode or COSEBI npz that holds vector and covariance together.

Copies arrays verbatim into two npz files, and for a `pure_eb_sweep` bundle does the
same to each member. `theta`, `theta_int`, `nmodes`, `version` and `*_scale_cut` go
to both files, so `reassemble` gives back the recorded arrays exactly.

    python tools/split_recorded_covariances.py pure_eb <pure_eb_data.npz> <out_data.npz> <out_covariance.npz>
    python tools/split_recorded_covariances.py cosebis <cosebis_modes_data.npz> <out_modes.npz> <out_covariance.npz>
    python tools/split_recorded_covariances.py pure_eb_sweep <pure_eb_sweep.tar.gz> <out_sweep.tar.gz> <out_covariance_sweep.tar.gz>
"""
import re
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np

PURE_SHARED = ("theta", "theta_int")
PURE_COVARIANCE = ("cov_pure_eb", "covariance_method", "n_samples")
COSEBI_SHARED = ("nmodes", "version", "fiducial_scale_cut", "full_scale_cut")
COSEBI_COVARIANCE = re.compile(r"(fiducial|full)_(cov|sigma_E|sigma_B)$")
MEMBER = re.compile(r"(?P<ver>.+)_(?P<blind>[A-C])_pure_eb_semianalytic\.npz$")


def split(recorded, statistic):
    """The (data, covariance) arrays of one recorded file, in its key order."""
    data, covariance = {}, {}
    for key in recorded.files:
        if statistic == "pure_eb":
            shared, is_covariance = key in PURE_SHARED, key in PURE_COVARIANCE
        else:
            shared, is_covariance = key in COSEBI_SHARED, bool(COSEBI_COVARIANCE.match(key))
        if shared or not is_covariance:
            data[key] = recorded[key]
        if shared or is_covariance:
            covariance[key] = recorded[key]
    return data, covariance


def reassemble(data, covariance):
    """The recorded arrays, from the two split files."""
    shared = set(data) & set(covariance)
    for key in shared:
        np.testing.assert_array_equal(data[key], covariance[key], key)
    return {**data, **covariance}


def split_file(statistic, recorded, out_data, out_covariance):
    with np.load(recorded) as source:
        data, covariance = split(source, statistic)
    np.savez(out_data, **data)
    np.savez(out_covariance, **covariance)


def split_sweep(bundle, out_data, out_covariance):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "in").mkdir()
        (tmp / "data").mkdir()
        (tmp / "covariance").mkdir()
        with tarfile.open(bundle) as tar:
            tar.extractall(tmp / "in", filter="data")
        for member in sorted((tmp / "in").glob("*.npz")):
            name = MEMBER.match(member.name)
            split_file("pure_eb", member,
                       tmp / "data" / f"{name['ver']}_pure_eb_data.npz",
                       tmp / "covariance" / f"{name['ver']}_{name['blind']}_pure_eb_covariance.npz")
        for out, side in ((out_data, "data"), (out_covariance, "covariance")):
            with tarfile.open(out, "w:gz") as tar:
                tar.add(tmp / side, arcname=".")


if __name__ == "__main__":
    kind, recorded, out_data, out_covariance = sys.argv[1:]
    if kind == "pure_eb_sweep":
        split_sweep(recorded, out_data, out_covariance)
    else:
        split_file(kind, recorded, out_data, out_covariance)
