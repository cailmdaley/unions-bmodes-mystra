"""Fiducial configuration-space scale cut.

A θ cut keeps the whole reporting bins that lie between its bounds, so no scale
outside the stated range is used. The fiducial cut is stated as 12′–83′. Those
bounds are edges 9 and 16 of the reporting grid ``geomspace(1', 250', 21)``,
11.997′ and 82.861′, quoted rounded to the arcminute, so the cut keeps
reporting bins 9–15 and spans exactly those edges. Functions take the reporting
grid as a ``decision_values.Grid``.
"""

import numpy as np

# The fiducial θ cut as stated, in arcmin.
STATED_CUT_ARCMIN = (12.0, 83.0)

# The stated bounds are bin edges rounded to the arcminute, so a kept bin's
# edges may miss them by that rounding (11.997' vs 12' is 2.5e-4 relative;
# 82.861' vs 83' is 1.7e-3). The nearest bin that is genuinely outside the cut
# misses by ~0.24 relative, so 1e-2 separates the two cases.
BIN_EDGE_RTOL = 1e-2


def reporting_grid_edges(grid):
    """Bin edges of the reporting grid, ``geomspace(1', 250', nbins + 1)``."""
    return np.geomspace(grid.min_sep, grid.max_sep, grid.nbins + 1)


def fiducial_bin_window(grid):
    """(first_bin, last_bin) of the fiducial cut: 0-based, inclusive.

    The reporting bins whose *edges* both lie between the stated bounds, which
    on the 20-bin grid are bins (9, 15).

    **Never resolve this window against bin mean angles.** The per-cut PTE
    grids ship a 20-element ``theta`` of bin *means*, and nearest-mean matching
    (``argmin(abs(theta - theta_min))``) gives the wrong window. Bin 8 spans
    9.105--11.997 arcmin with a mean of 10.614, bin 9 spans 11.997--15.807 with
    a mean of 13.988, so 11.997 snaps to bin 8: the window becomes bins (8, 15),
    i.e. 9.1--83 arcmin, one bin wider than the declared cut, with nu = 2*8 = 16
    rather than 14 in the joint xi_tot^B test.
    """
    theta_min, theta_max = STATED_CUT_ARCMIN
    edges = reporting_grid_edges(grid)
    inside = (edges[:-1] >= theta_min * (1.0 - BIN_EDGE_RTOL)) & (
        edges[1:] <= theta_max * (1.0 + BIN_EDGE_RTOL)
    )
    bins = np.flatnonzero(inside)
    return int(bins[0]), int(bins[-1])


def fiducial_theta_cut(grid):
    """(theta_min, theta_max) of the fiducial cut, in arcmin: the outer edges of
    the kept reporting bins, (11.997, 82.861) on the 20-bin grid."""
    first, last = fiducial_bin_window(grid)
    edges = reporting_grid_edges(grid)
    return float(edges[first]), float(edges[last + 1])
