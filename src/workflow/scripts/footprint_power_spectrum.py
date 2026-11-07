"""Normalised power spectrum of the survey footprint, from the comprehensive catalog.

The footprint is the HEALPix nside-4096 map of pixels holding at least one object
that passes the spatially structured cuts of the comprehensive catalog; per-galaxy
quality cuts play no part. It is the footprint the size-cut, initial and
relaxed-flags catalogs share, 14 123 720 pixels. The spectrum is normalised so
that Σ_ℓ (2ℓ+1)/4π C_ℓ, in deg², equals the footprint area, and written as
`ell C_ell` rows for CosmoCov.

    python footprint_power_spectrum.py \
        --comprehensive-catalog <comprehensive.hdf5> --out-file <spectrum.txt>
"""

import argparse

import h5py
import healpy as hp
import numpy as np

NSIDE = 4096
CHUNK = 500_000


def in_footprint(data, flags):
    """The spatially structured cuts, on the `data` and `data_ext` groups."""
    return (
        (data["overlap"] == 1)
        & (data["IMAFLAGS_ISO"] == 0)
        & (data["N_EPOCH"] >= 2)
        & (flags["4_Stars"] == 0)
        & (flags["8_Manual"] == 0)
        & (flags["64_r"] == 0)
        & (flags["1024_Maximask"] == 0)
        & (flags["npoint3"] >= 3)
    )


def footprint_map(path, nside=NSIDE):
    """Binary HEALPix map of the pixels occupied after the cuts."""
    footprint = np.zeros(hp.nside2npix(nside), dtype=np.uint8)
    with h5py.File(path, "r") as catalog:
        for start in range(0, catalog["data"].shape[0], CHUNK):
            data = catalog["data"][start : start + CHUNK]
            flags = catalog["data_ext"][start : start + CHUNK]
            keep = in_footprint(data, flags)
            footprint[hp.ang2pix(nside, data["RA"][keep], data["Dec"][keep], lonlat=True)] = 1
    return footprint


def normalised_spectrum(footprint, nside=NSIDE):
    """The footprint's C_ℓ up to 3 nside − 1, with its area as the normalisation."""
    area = footprint.sum() * hp.nside2pixarea(nside, degrees=True)
    cl = hp.anafast(footprint, lmax=3 * nside - 1)
    ell = np.arange(len(cl))
    integral = np.sum((2 * ell + 1) / (4 * np.pi) * cl) / (np.pi / 180) ** 2
    return ell, cl * area / integral, area


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--comprehensive-catalog", required=True)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    footprint = footprint_map(a.comprehensive_catalog)
    ell, spectrum, area = normalised_spectrum(footprint.astype(float))
    np.savetxt(a.out_file, np.column_stack((ell, spectrum)), fmt=["%d", "%.10e"])
    print(f"{int(footprint.sum())} pixels, {area:.2f} deg2 -> {a.out_file}")


if __name__ == "__main__":
    main()
