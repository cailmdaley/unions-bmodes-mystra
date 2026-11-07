"""Audit the sky the data and one mock catalogue occupy, at HEALPix nside 4096.

For each catalogue, its row count, occupied pixels and their area, then the two
footprints' overlap and each one's exclusive area, as text lines:

    data: N=<rows>, occupied pix=<pixels>, area=<deg2> deg2
    mock: N=<rows>, occupied pix=<pixels>, area=<deg2> deg2
    overlap = <deg2> deg2
    data only = <deg2> deg2
    mock only = <deg2> deg2

    python mock_footprint_audit.py --cat-config <cat_config.yaml> \
        --data-version SP_v1.4.6.3_leak_corr \
        --mock-version SP_v1.4.6.3_glass_mock_v2 --out-file <audit.txt>
"""

import argparse
import os

import healpy as hp
import numpy as np
import yaml
from astropy.io import fits

NSIDE = 4096


def occupied_pixels(entry):
    """The row count and sorted unique pixels of a cat_config entry's shear catalog."""
    path = os.path.join(entry["subdir"], entry["shear"]["path"])
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        return len(data), np.unique(hp.ang2pix(NSIDE, data["ra"], data["dec"], lonlat=True))


def audit_lines(data_rows, data_pix, mock_rows, mock_pix):
    pixel_area = hp.nside2pixarea(NSIDE, degrees=True)
    overlap = np.intersect1d(data_pix, mock_pix, assume_unique=True).size
    return [
        f"data: N={data_rows}, occupied pix={data_pix.size}, area={data_pix.size * pixel_area:.2f} deg2",
        f"mock: N={mock_rows}, occupied pix={mock_pix.size}, area={mock_pix.size * pixel_area:.2f} deg2",
        f"overlap = {overlap * pixel_area:.2f} deg2",
        f"data only = {(data_pix.size - overlap) * pixel_area:.2f} deg2",
        f"mock only = {(mock_pix.size - overlap) * pixel_area:.2f} deg2",
    ]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cat-config", required=True)
    ap.add_argument("--data-version", required=True)
    ap.add_argument("--mock-version", required=True)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    with open(a.cat_config) as stream:
        entries = yaml.safe_load(stream)
    data_rows, data_pix = occupied_pixels(entries[a.data_version.removesuffix("_leak_corr")])
    mock_rows, mock_pix = occupied_pixels(entries[a.mock_version])
    lines = audit_lines(data_rows, data_pix, mock_rows, mock_pix)
    with open(a.out_file, "w") as stream:
        stream.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
