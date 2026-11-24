"""High-resolution TreeCorr ξ± of one catalog version on one node.

Bins this fine shrink TreeCorr's tree-opening criterion, so the pair traversal
is long: about 720 core-hours for 10,000 bins. The script reads the
catalog entry straight from cat_config and takes its thread count from
OMP_NUM_THREADS, which `lc` sets from the recipe's cpus.

    python run_2pcf_highres.py \
        --cat-config src/cosmo_val/cat_config.yaml --version SP_v1.4.6.3_leak_corr \
        --nbins 10000 --min-sep 0.5 --max-sep 300.0 --out-file <xi.txt>
"""
import argparse
import copy
import os
import time

import numpy as np
import treecorr
import yaml
from astropy.io import fits
# Shear response, 1 for the ShapePipe catalogs.
R = 1.0
LEAK_SUFFIX = "_leak_corr"


def log(message):
    """Print with a timestamp."""
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def shear_config(cat_config, version):
    """The catalog path and columns of a version, with its `_leak_corr` columns."""
    with open(cat_config) as stream:
        entries = yaml.safe_load(stream)
    entry = copy.deepcopy(entries[version.removesuffix(LEAK_SUFFIX)]["shear"])
    if version.endswith(LEAK_SUFFIX):
        entry["e1_col"] = entry["e1_col_corrected"]
        entry["e2_col"] = entry["e2_col_corrected"]
    subdir = entries[version.removesuffix(LEAK_SUFFIX)]["subdir"]
    entry["path"] = os.path.join(subdir, entry["path"])
    return entry


def load_catalog(shear):
    """Load the shear catalog and subtract the weighted mean ellipticity."""
    with fits.open(shear["path"], memmap=True) as hdul:
        data = hdul[1].data
        ra = np.array(data["ra"], dtype=np.float64)
        dec = np.array(data["dec"], dtype=np.float64)
        e1 = np.array(data[shear["e1_col"]], dtype=np.float64)
        e2 = np.array(data[shear["e2_col"]], dtype=np.float64)
        w = np.array(data[shear["w_col"]], dtype=np.float64)
    log(f"{len(ra):,} galaxies loaded from {shear['path']}")
    g1 = (e1 - np.average(e1 / R, weights=w)) / R
    g2 = (e2 - np.average(e2 / R, weights=w)) / R
    return ra, dec, g1, g2, w


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cat-config", required=True)
    ap.add_argument("--version", required=True, help="Catalog version key in cat_config")
    ap.add_argument("--nbins", type=int, required=True)
    ap.add_argument("--min-sep", type=float, required=True, help="arcmin")
    ap.add_argument("--max-sep", type=float, required=True, help="arcmin")
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)

    start = time.time()
    threads = int(os.environ["OMP_NUM_THREADS"])
    log(f"{a.version}: {a.nbins:,} bins, [{a.min_sep}, {a.max_sep}] arcmin, "
        f"{threads} threads")

    ra, dec, g1, g2, w = load_catalog(shear_config(a.cat_config, a.version))
    cat = treecorr.Catalog(
        ra=ra, dec=dec, g1=g1, g2=g2, w=w, ra_units="degrees", dec_units="degrees",
    )
    del ra, dec, g1, g2, w

    gg = treecorr.GGCorrelation(
        min_sep=a.min_sep, max_sep=a.max_sep, nbins=a.nbins,
        sep_units="arcminutes", verbose=2,
    )
    gg.process(cat, num_threads=threads)
    gg.write(a.out_file)
    log(f"Wrote {a.out_file} after {(time.time() - start) / 3600:.1f} h")


if __name__ == "__main__":
    main()
