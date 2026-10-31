"""Fiducial angular scale cut, resolved under the scale_cut_convention decision.

The paper's fiducial configuration-space cut exists in two forms that disagree
at the third digit. The published pipeline cut at *edges* of the 21-bin
reporting grid ``geomspace(1', 250', 21)`` — indices 9 and 16, i.e.
(11.997…, 82.861…) — while the paper text and ``config.yaml`` round those to
[12, 83]. ``scale_cut_convention`` in the ASTRA spec selects between them; this
module is the single place the choice is resolved.

Applies to θ cuts only: the harmonic ℓ cut (300–1600) has no grid/rounding
ambiguity.
"""

import numpy as np

CONVENTIONS = ("grid_bin_edges", "rounded_arcmin")

# Indices into the reporting grid's bin edges that the published pipeline used
# for the fiducial cut.
FIDUCIAL_EDGE_INDICES = (9, 16)


def reporting_grid_edges(config):
    """Bin edges of the reporting grid, ``geomspace(1', 250', nbins + 1)``."""
    fid = config["fiducial"]
    return np.geomspace(
        float(fid["min_sep"]), float(fid["max_sep"]), int(fid["nbins"]) + 1
    )


def resolve_fiducial_scale_cuts(config, convention, config_key=None):
    """(theta_min, theta_max) of the fiducial cut, in arcmin.

    Parameters
    ----------
    config : dict
        Parsed bmodes ``config.yaml``.
    convention : {'grid_bin_edges', 'rounded_arcmin'}
        See module docstring.
    config_key : str, optional
        Under ``rounded_arcmin``, read the pair from this key of the
        ``fiducial`` block (e.g. ``fiducial_xip_scale_cut``) instead of the
        default ``fiducial_min_scale`` / ``fiducial_max_scale`` scalars.
        Ignored under ``grid_bin_edges``, which is a property of the grid.
    """
    if convention == "grid_bin_edges":
        edges = reporting_grid_edges(config)
        lo, hi = (edges[i] for i in FIDUCIAL_EDGE_INDICES)
    elif convention == "rounded_arcmin":
        fid = config["fiducial"]
        if config_key is None:
            lo, hi = fid["fiducial_min_scale"], fid["fiducial_max_scale"]
        else:
            lo, hi = fid[config_key]
    else:
        raise ValueError(
            f"unknown scale-cut convention {convention!r}; expected one of {CONVENTIONS}"
        )
    return float(lo), float(hi)


def add_scale_cut_convention_argument(parser):
    """Add the required ``--scale-cut-convention`` flag.

    Deliberately has no default: the ASTRA spec is the single source of truth
    for which convention a given universe adopts.
    """
    parser.add_argument(
        "--scale-cut-convention",
        required=True,
        choices=CONVENTIONS,
        help=(
            "How the fiducial theta cut is defined: 'grid_bin_edges' uses edges "
            "9/16 of the geomspace(1, 250, 21) reporting grid (what the published "
            "pipeline computed); 'rounded_arcmin' uses the rounded config values "
            "([12, 83]) as printed in the paper text."
        ),
    )


# A bin's edges may miss the resolved cut by the rounding between the two
# conventions (11.997' vs 12' is 2.5e-4 relative; 82.861' vs 83' is 1.7e-3).
# The nearest bin that is genuinely outside the window misses by ~0.24
# relative, so 1e-2 separates the two cases with three decades to spare.
BIN_EDGE_RTOL = 1e-2


def resolve_fiducial_bin_window(config, convention, config_key=None):
    """(first_bin, last_bin) of the fiducial cut: 0-based, inclusive.

    The reporting bins whose *edges* both lie inside the resolved
    (theta_min, theta_max), against ``reporting_grid_edges(config)``. Under
    either convention the 12--83 arcmin cut resolves to bins (9, 15) — seven
    bins, spanning edges 9 through 16 of ``geomspace(1', 250', 21)``.

    **Never resolve this window against bin mean angles.** The per-cut PTE
    grids ship a 20-element ``theta`` of bin *means*, and the code this
    function replaces did ``argmin(abs(theta - theta_min))`` against it. Bin 8
    spans 9.105--11.997 arcmin with a mean of 10.614, bin 9 spans
    11.997--15.807 with a mean of 13.988, so 11.997 snapped to bin 8: the
    window silently became bins (8, 15), i.e. 9.1--83 arcmin, one bin wider
    than the cut the paper declares. Every published and revised fiducial
    pure-mode PTE (xi+^B, xi-^B, xi_tot^B) was computed on that wider window,
    and the dof-penalised combined PTE inherited nu = 2*8 = 16 instead of 14.
    """
    theta_min, theta_max = resolve_fiducial_scale_cuts(config, convention, config_key)
    edges = reporting_grid_edges(config)
    left, right = edges[:-1], edges[1:]
    inside = (left >= theta_min * (1.0 - BIN_EDGE_RTOL)) & (
        right <= theta_max * (1.0 + BIN_EDGE_RTOL)
    )
    bins = np.flatnonzero(inside)
    if bins.size == 0:
        raise ValueError(
            f"no reporting bin lies inside ({theta_min}, {theta_max}) arcmin"
        )
    return int(bins[0]), int(bins[-1])
