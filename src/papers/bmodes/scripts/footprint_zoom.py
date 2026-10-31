"""Paper II Fig. 1 for a catalog comparison: galaxy positions around MegaCam exposure 2278619.

Two stacked 2:1 panels at single-column width, each 144' x 72' centred on the exposure's
mosaic centre (248.96, +54.50), east left, no axes. Top: the reference catalog with the 40 CCD
outlines; bottom: the fiducial catalog. Dots are galaxy positions from the cut catalogs.

    python footprint_zoom.py --top <reference.fits> --bottom <fiducial.fits> \
        --ccds <exp_2278619_ccds.json> --out-file <figure.pdf>
"""
import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.collections
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from plotting_utils import PAPER_MPLSTYLE

W1 = 8.8 / 2.54          # column width [in]
HALF_H = 36.0            # panel half-height [arcmin]; the mosaic spans +-30'
HALF_W = 2 * HALF_H
RA_C, DEC_C = 248.9598, 54.5048
BLUE, RED = "#2a78d6", "crimson"
DOT_PT, DOT_ALPHA = 1.5, 0.3
BAR = 10.0               # scale bar [arcmin]

w = WCS(naxis=2)
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w.wcs.crval = [RA_C, DEC_C]
w.wcs.crpix = [1, 1]
w.wcs.cdelt = [-1 / 60, 1 / 60]  # one pixel = 1', x grows westward so east is left


def off(ra, dec):
    return w.world_to_pixel_values(ra, dec)


from astropy.io import fits


def ccd_segments(polys, tol=0.5):
    """CCD edges, with each inter-chip gap drawn once along its centre line."""
    edges = [(i, np.array([p[k], p[(k + 1) % 4]])) for i, p in enumerate(polys) for k in range(4)]
    used, segs = set(), []
    for a, (i, ea) in enumerate(edges):
        if a in used:
            continue
        for b in range(a + 1, len(edges)):
            j, eb = edges[b]
            if b in used or j == i or np.hypot(*(ea.mean(0) - eb.mean(0))) > tol:
                continue
            if np.hypot(*(ea[0] - eb[0])) > np.hypot(*(ea[0] - eb[1])):
                eb = eb[::-1]
            ea = 0.5 * (ea + eb)
            used.add(b)
            break
        used.add(a)
        segs.append(ea)
    return segs


def frame(ax):
    ax.set_xlim(-HALF_W, HALF_W)
    ax.set_ylim(-HALF_H, HALF_H)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.6)


def scale_bar(ax):
    """Solid 10' bar with its label centred above, together in a small white box."""
    sb = AnchoredSizeBar(
        ax.transData, BAR, rf"${BAR:g}'$", loc="lower left", pad=0.25, borderpad=0.5, sep=2,
        frameon=True, size_vertical=0, label_top=True, fontproperties=FontProperties(size=LABEL),
    )
    sb.size_bar.get_children()[0].set_linewidth(1.2)
    sb.patch.set_facecolor("white")
    sb.patch.set_edgecolor("none")
    sb.zorder = 7
    ax.add_artist(sb)


def positions(path):
    """Galaxy offsets [arcmin] inside the panel."""
    with fits.open(path, memmap=True) as h:
        ra, dec = (np.asarray(h[1].data[c], dtype=np.float64) for c in ("RA", "Dec"))
    box = (np.abs(dec - DEC_C) < 0.8) & (np.abs(ra - RA_C) < 2.4)
    x, y = off(ra[box], dec[box])
    s = (np.abs(x) < HALF_W + 1) & (np.abs(y) < HALF_H + 1)
    return x[s], y[s]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--top", required=True)
    ap.add_argument("--bottom", required=True)
    ap.add_argument("--ccds", required=True)
    ap.add_argument("--out-file", required=True)
    a = ap.parse_args(argv)
    plt.style.use(PAPER_MPLSTYLE)
    plt.rcParams.update({"figure.autolayout": False})
    global LABEL
    LABEL = plt.rcParams["axes.labelsize"]
    ex = json.load(open(a.ccds))
    ccds = [np.column_stack(off(np.array(c["ra"]), np.array(c["dec"]))) for c in ex["ccds"]]
    fig, axes = plt.subplots(2, 1, figsize=(W1, W1 * 0.5 * 2 + 0.06), gridspec_kw=dict(hspace=0.04))
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    for ax, path in zip(axes, (a.top, a.bottom)):
        x, y = positions(path)
        ax.scatter(x, y, s=DOT_PT**2, lw=0, color=BLUE, alpha=DOT_ALPHA, rasterized=True)
        frame(ax)
    axes[0].add_collection(matplotlib.collections.LineCollection(
        ccd_segments(ccds), colors=RED, linewidths=0.8, capstyle="projecting", zorder=4))
    scale_bar(axes[0])
    fig.savefig(a.out_file, dpi=600)


if __name__ == "__main__":
    main()
