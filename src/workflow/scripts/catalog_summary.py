"""What the covariances and the GLASS mocks need from a shear catalog.

Reads each catalog once and writes one npz holding, per catalog key K:

- ``K/pix1024``, ``K/wsum1024``: the occupied nside-1024 pixels (RING) and the
  summed weight in each, the weighted count map the NaMaster Cℓ covariance takes
  as its window;
- ``K/occupied4096``: the nside-4096 pixels holding at least one galaxy, as
  ``np.packbits`` of the full-sky boolean map, the GLASS footprint;
- ``K/n``, ``K/sum_w``, ``K/sum_w2``: galaxy count, Σw, Σw²;
- ``K/sum_w2e2``: Σw²(e1² + e2²)/2 on the raw ellipticities, and
  ``K/sum_w2e2_leak_corr`` on the leakage-corrected ones when the catalog has
  them.

The survey numbers follow: the weight-effective count (Σw)²/Σw², the
per-component σ_e² = sum_w2e2 / Σw², and the footprint area from the nside-4096
occupancy (:func:`survey`).

    python catalog_summary.py --catalog SP_v1.4.6.3=cat.fits [...] --out-file summary.npz
"""

import argparse

import fitsio
import healpy as hp
import numpy as np

NSIDE_WINDOW = 1024
NSIDE_FOOTPRINT = 4096


def summarize(path):
    columns = fitsio.FITS(path)[1].get_colnames()
    weight = "w_des" if "w_des" in columns else "w"
    names = ["RA", "Dec", weight, "e1", "e2"]
    corrected = "e1_leak_corrected" in columns
    if corrected:
        names += ["e1_leak_corrected", "e2_leak_corrected"]
    d = fitsio.read(path, columns=names)
    ra, dec = d["RA"].astype(np.float64), d["Dec"].astype(np.float64)
    w = d[weight].astype(np.float64)

    # Pixelised exactly as sp_validation.pseudo_cl.get_n_gal_map does.
    pix = hp.ang2pix(NSIDE_WINDOW, (90.0 - dec) * np.pi / 180.0, ra * np.pi / 180.0)
    pix1024, inverse = np.unique(pix, return_inverse=True)
    wsum1024 = np.bincount(inverse, weights=w)

    occupied = np.zeros(hp.nside2npix(NSIDE_FOOTPRINT), dtype=bool)
    occupied[hp.ang2pix(NSIDE_FOOTPRINT, ra, dec, lonlat=True)] = True

    w2 = w**2

    def sum_w2e2(e1, e2):
        e1, e2 = e1.astype(np.float64), e2.astype(np.float64)
        return np.sum(w2 * (e1**2 + e2**2) / 2)

    out = {
        "pix1024": pix1024.astype(np.int64),
        "wsum1024": wsum1024,
        "occupied4096": np.packbits(occupied),
        "n": np.int64(len(w)),
        "sum_w": np.sum(w),
        "sum_w2": np.sum(w2),
        "sum_w2e2": sum_w2e2(d["e1"], d["e2"]),
    }
    if corrected:
        out["sum_w2e2_leak_corr"] = sum_w2e2(d["e1_leak_corrected"], d["e2_leak_corrected"])
    return out


def load(path, key):
    """One catalog's entries from a summary, without the ``K/`` prefix."""
    with np.load(path) as z:
        return {name.split("/", 1)[1]: z[name] for name in z.files if name.split("/", 1)[0] == key}


def window(entry):
    """The full-sky nside-1024 weighted count map."""
    n_gal = np.zeros(hp.nside2npix(NSIDE_WINDOW))
    n_gal[entry["pix1024"]] = entry["wsum1024"]
    return n_gal


def footprint(entry):
    """The full-sky nside-4096 occupancy map, as float64 0/1."""
    bits = np.unpackbits(entry["occupied4096"], count=hp.nside2npix(NSIDE_FOOTPRINT))
    return bits.astype(np.float64)


def survey(entry, leak_corr=False):
    """Footprint area (deg²), weight-effective count, per-component σ_e."""
    occupied = int(np.unpackbits(entry["occupied4096"]).sum())
    area = float(occupied * hp.nside2pixarea(NSIDE_FOOTPRINT, degrees=True))
    n_eff = float(entry["sum_w"] ** 2 / entry["sum_w2"])
    s = entry["sum_w2e2_leak_corr" if leak_corr else "sum_w2e2"]
    return area, n_eff, float(np.sqrt(s / entry["sum_w2"]))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog", action="append", required=True, metavar="KEY=PATH")
    p.add_argument("--out-file", required=True)
    a = p.parse_args(argv)
    arrays = {}
    for spec in a.catalog:
        key, path = spec.split("=", 1)
        entry = summarize(path)
        arrays.update({f"{key}/{name}": value for name, value in entry.items()})
        area, n_eff, sigma_e = survey(entry)
        print(f"{key}: {entry['n']} galaxies, area {area:.2f} deg², n_eff {n_eff:.0f}, σ_e {sigma_e:.6f}")
    with open(a.out_file, "wb") as f:
        np.savez_compressed(f, **arrays)


if __name__ == "__main__":
    main()
