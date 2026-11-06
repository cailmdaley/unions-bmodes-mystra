"""Shared plotting utilities for claims scripts."""

from pathlib import Path

import matplotlib.scale as mscale
import matplotlib.ticker as ticker
import matplotlib.transforms as mtransforms
import numpy as np
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle

# Shared constants
# Resolve the paper style relative to this module (.resolve() follows the
# papers/bmodes/scripts symlink to the real workflow/scripts location, so both
# import paths land on workflow/matplotlib_config/); survives project moves.
_MPL_CONFIG = Path(__file__).resolve().parent.parent / "matplotlib_config"
PAPER_MPLSTYLE = str(_MPL_CONFIG / "paper.mplstyle")

# Colour encodes the *estimator family* in every B-mode data-vector figure
# (paper Figs. 2, 3, 4, 6): each family has one hue, and a second hue appears
# only for a peer series that shares the panel.  Within a family the filled
# marker is the fiducial series and the open marker its peer, so the
# encoding survives greyscale.
COLOR_COSEBIS = "#7b3fa0"           # B_n at the fiducial cuts (filled)
COLOR_COSEBIS_ALL = "#c9a6dd"       # B_n over all scales (open)
COLOR_XI_B = "crimson"              # xi_+^B (filled); xi^B in Fig. 2
COLOR_XI_B_MINUS = "#e8792b"        # xi_-^B (open), where it shares a panel
COLOR_CL_BB = "#008080"             # C_l^BB (filled)
COLOR_CL_EB = "#4d4d4d"             # C_l^EB (open)
# Fig. 2 decomposition: the modes of xi_pm, B in the family hue.
COLOR_TOTAL = "k"
COLOR_E = "#2a78d6"
COLOR_B = COLOR_XI_B
COLOR_AMB = "#9a9a9a"

# One mark for every B-mode data-vector figure (paper Figs. 2, 3, 4, 6):
# small points, thin bars, subtle caps.  Series on the same axes differ by
# hue (see the family palette above) and, redundantly, by marker fill.
DATA_VECTOR_MARK = {
    "ls": "none",
    "markersize": 3.2,
    "markeredgewidth": 0.8,
    "elinewidth": 0.7,
    "capsize": 1.5,
    "capthick": 0.7,
}

# Scales outside the fiducial selection, shaded the same way in every figure.
EXCLUDED_SHADE = {"color": "gray", "alpha": 0.1, "lw": 0, "zorder": 0}
FIG_WIDTH_FULL = 7.24  # Two-column A&A format (textwidth)
FIG_WIDTH_SINGLE = 3.54  # Single-column A&A format (columnwidth)
MARKER_STYLES = ["o", "s", "D", "^"]

# Categorical series colours, slots 1-2 of the data-viz reference palette
# (blue, orange).  Chosen over a generated husl pair because they clear the
# all-pairs colour-vision and normal-vision separation floors with room to
# spare (CVD ΔE 24.7, normal ΔE 33.6 against the light surface) and differ in
# lightness, so they survive greyscale printing; MARKER_STYLES carries the
# secondary encoding on top.
CATEGORICAL_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]

# Scatter-with-uncertainty marks: no connecting line (these are independent
# per-bin estimates, not a series in a continuous variable), thin bars, a
# surface-coloured ring so overlapping series stay separable.
SCATTER_ERRORBAR = {
    "ls": "none",
    "ms": 4,
    "elinewidth": 0.8,
    "capsize": 0,
    "markeredgewidth": 0.5,
    "markeredgecolor": "white",
}

# Default errorbar styling for version comparison plots
ERRORBAR_DEFAULTS = {
    "markeredgewidth": 0.5,
    "markersize": 4,
    "capsize": 2,
    "capthick": 0.8,
    "elinewidth": 0.8,
}

# Version-comparison figures.
X_OFFSETS = [-0.132, -0.044, 0.044, 0.132]
VERSION_ALPHA = {"fiducial": 1.0, "comparison": 1.0}
VERSION_BOX = {
    "edge_color": "0.3",
    "edge_linewidth": 0.7,
    "fiducial_line_color": "black",
    "fiducial_line_width": 1.0,
}


PTE_GRADIENT_LEVELS = 4096
def PTE_COLORBAR_TICKS(threshold):
    """Return the threshold-focused PTE colorbar tick locations."""
    return [threshold, 0.5, 1 - threshold]


def PTE_COLORBAR_TICKLABELS(threshold):
    """Return labels for the threshold-focused PTE colorbar ticks."""
    return [f"{value:g}" for value in PTE_COLORBAR_TICKS(threshold)]


def make_pte_colormap(
    threshold, gradient_range=(0.15, 0.85), n_gradient=PTE_GRADIENT_LEVELS
):
    """Create a discrete PTE colormap with exact threshold breaks.

    Values below ``threshold`` and above ``1 - threshold`` remain solid
    colors. The interval between them uses a high-resolution gradient so the
    threshold boundary is precise without losing the discrete fail/pass cue.
    """
    vlag = sns.color_palette("vlag", as_cmap=True)
    solid_blue = vlag(0.0)
    solid_red = vlag(1.0)
    g_lo, g_hi = gradient_range
    gradient_colors = [vlag(val) for val in np.linspace(g_lo, g_hi, n_gradient)]
    cmap = ListedColormap(
        [solid_blue] + gradient_colors + [solid_red],
        name="pte_threshold",
    )
    cmap.set_bad(color="lightgray")
    return cmap


def make_pte_norm(threshold, n_gradient=PTE_GRADIENT_LEVELS):
    """Create exact threshold normalization for PTE heatmaps."""
    low = threshold
    high = 1 - threshold
    boundaries = np.concatenate(
        (
            [0.0],
            np.linspace(low, high, n_gradient + 1),
            [1.0],
        )
    )
    return BoundaryNorm(boundaries, ncolors=n_gradient + 2, clip=True)


def format_pte_colorbar(cbar, threshold):
    """Apply consistent threshold-focused styling to PTE colorbars."""
    cbar.set_ticks(PTE_COLORBAR_TICKS(threshold))
    cbar.set_ticklabels(PTE_COLORBAR_TICKLABELS(threshold))
    cbar.minorticks_off()
    if cbar.solids is not None:
        cbar.solids.set_edgecolor("face")
        cbar.solids.set_linewidth(0.0)


class SquareRootScale(mscale.ScaleBase):
    """Square root scale for x-axis (matches bandpower binning)."""

    name = "squareroot"

    def __init__(self, axis, **kwargs):
        mscale.ScaleBase.__init__(self, axis, **kwargs)

    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(ticker.AutoLocator())
        axis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        axis.set_minor_locator(ticker.NullLocator())
        axis.set_minor_formatter(ticker.NullFormatter())

    def limit_range_for_scale(self, vmin, vmax, minpos):
        return max(0.0, vmin), vmax

    class SquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def transform_non_affine(self, a):
            return np.array(a) ** 0.5

        def inverted(self):
            return SquareRootScale.InvertedSquareRootTransform()

    class InvertedSquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def transform_non_affine(self, a):
            return np.array(a) ** 2

        def inverted(self):
            return SquareRootScale.SquareRootTransform()

    def get_transform(self):
        return self.SquareRootTransform()


# Register at import time
mscale.register_scale(SquareRootScale)


def version_label(version, version_labels):
    """Get the human-readable label for a version.

    Parameters
    ----------
    version : str
        Version string (e.g., "SP_v1.4.11.2_leak_corr").
    version_labels : dict
        Mapping from version strings to human-readable labels.

    Returns
    -------
    str
        Human-readable label, or cleaned version string if not in mapping.
    """
    return version_labels.get(
        version, version.replace("SP_", "").replace("_leak_corr", "")
    )


def find_fiducial_index(datasets, fiducial_version, key="version"):
    """Find index of fiducial version in datasets list.

    Parameters
    ----------
    datasets : list of dict
        List of dataset dictionaries with version key.
    fiducial_version : str
        Version string to find.
    key : str
        Key to use for matching (default: "version").

    Returns
    -------
    int
        Index of fiducial version, or 0 if not found.
    """
    return next((i for i, d in enumerate(datasets) if d[key] == fiducial_version), 0)


def get_version_alpha(version, fiducial_version):
    """Get alpha value for version based on whether it's fiducial.

    Parameters
    ----------
    version : str
        Version string.
    fiducial_version : str
        Fiducial version string for comparison.

    Returns
    -------
    float
        Alpha value from ``VERSION_ALPHA``.
    """
    return VERSION_ALPHA["fiducial" if version == fiducial_version else "comparison"]


def draw_normalized_version_box(ax, x_left, x_right, y_vals, fiducial_val, style):
    """Draw a single version box for normalized data.

    Draws a rectangle spanning the range of values +/- 1 sigma (where sigma=1
    by construction for normalized data), with a fiducial reference line.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    x_left, x_right : float
        Left and right edges of box.
    y_vals : array-like
        Normalized y values for all versions at this bin.
    fiducial_val : float
        Fiducial version's y value for the reference line.
    style : dict
        Styling dict, ``VERSION_BOX``.
    """
    y_vals = np.asarray(y_vals)
    box_bottom = y_vals.min()
    box_top = y_vals.max()

    rect = Rectangle(
        (x_left, box_bottom),
        x_right - x_left,
        box_top - box_bottom,
        facecolor="none",
        edgecolor=style["edge_color"],
        linewidth=style["edge_linewidth"],
        zorder=1,
    )
    ax.add_patch(rect)

    ax.hlines(
        fiducial_val,
        x_left,
        x_right,
        colors=style["fiducial_line_color"],
        linewidth=style["fiducial_line_width"],
        zorder=1,
    )


def draw_normalized_boxes_log_scale(
    ax, x_centers, datasets, y_norm_key, fiducial_idx, x_offset_range
):
    """Draw version boxes for log-scale x-axis (e.g., theta).

    For each bin, draws a box spanning the range of values across versions,
    with a horizontal line at the fiducial version's value.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    x_centers : array-like
        Bin centers (theta values).
    datasets : list of dict
        Dataset dicts containing y_norm_key.
    y_norm_key : str
        Key for normalized y values in datasets.
    fiducial_idx : int
        Index of fiducial version in datasets.
    x_offset_range : tuple
        (min, max) multiplicative offsets for data points. Box covers
        this range plus 10% padding.
    """
    style = VERSION_BOX
    offset_min, offset_max = x_offset_range
    padding = 0.25 * (offset_max - offset_min)
    box_left_factor = offset_min - padding
    box_right_factor = offset_max + padding

    for i, x_i in enumerate(x_centers):
        y_vals = np.array([data[y_norm_key][i] for data in datasets])
        x_left = x_i * box_left_factor
        x_right = x_i * box_right_factor
        draw_normalized_version_box(
            ax, x_left, x_right, y_vals, y_vals[fiducial_idx], style
        )


def draw_normalized_boxes_linear_scale(
    ax, x_centers, datasets, y_norm_key, fiducial_idx, x_offsets
):
    """Draw version boxes for linear-scale x-axis (e.g., COSEBIS modes).

    For each bin, draws a box spanning the range of values across versions,
    with a horizontal line at the fiducial version's value.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    x_centers : array-like
        Bin centers (mode numbers).
    datasets : list of dict
        Dataset dicts containing y_norm_key.
    y_norm_key : str
        Key for normalized y values in datasets.
    fiducial_idx : int
        Index of fiducial version in datasets.
    x_offsets : array-like
        Array of x-offsets used for jittering data points. Box width covers
        this range plus 25% padding.
    """
    style = VERSION_BOX
    x_offsets = np.asarray(x_offsets)
    box_half_width = np.max(np.abs(x_offsets)) * 1.7

    for i, x_i in enumerate(x_centers):
        y_vals = np.array([data[y_norm_key][i] for data in datasets])
        x_left = x_i - box_half_width
        x_right = x_i + box_half_width
        draw_normalized_version_box(
            ax, x_left, x_right, y_vals, y_vals[fiducial_idx], style
        )


def draw_normalized_boxes_ell_scale(
    ax,
    ell,
    ell_widths,
    datasets,
    y_norm_key,
    fiducial_idx,
    jitter_fraction,
    n_versions,
):
    """Draw version boxes for ell-space plots with variable bin widths.

    For each multipole bin, draws a box spanning the range of values across
    versions, with a horizontal line at the fiducial version's value.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    ell : array-like
        Multipole bin centers.
    ell_widths : array-like
        Width of each ell bin.
    datasets : list of dict
        Dataset dicts containing y_norm_key.
    y_norm_key : str
        Key for normalized y values in datasets.
    fiducial_idx : int
        Index of fiducial version in datasets.
    jitter_fraction : float
        Jitter fraction used for offsetting data points.
    n_versions : int
        Number of versions being plotted.
    """
    style = VERSION_BOX
    # Max jitter is ((n-1)/2) * jitter_fraction; add 15% padding
    box_half_width_factor = ((n_versions - 1) / 2) * jitter_fraction * 1.4

    for i, ell_i in enumerate(ell):
        y_vals = np.array([data[y_norm_key][i] for data in datasets])
        half_width = ell_widths[i] * box_half_width_factor
        x_left = ell_i - half_width
        x_right = ell_i + half_width
        draw_normalized_version_box(
            ax, x_left, x_right, y_vals, y_vals[fiducial_idx], style
        )


def extract_version_number(version_string):
    """Extract short version number from full version string.

    Examples:
        SP_v1.4.5_leak_corr -> v1.4.5
        SP_v1.4.11.2_leak_corr -> v1.4.11.2
        SP_v1.4.5 -> v1.4.5

    Parameters
    ----------
    version_string : str
        Full version string (e.g., "SP_v1.4.5_leak_corr").

    Returns
    -------
    str
        Short version number (e.g., "v1.4.5").
    """
    import re

    match = re.search(r"(v[\d.]+)", version_string)
    return match.group(1) if match else version_string


def iter_version_figures(version_labels, fiducial_version):
    """Iterate over the 9 figure specs for per-version data vector plots.

    Produces specs for:
    - 1 paper figure: fiducial version, leak-corrected, no title
    - 8 dashboard figures: all versions × (corrected + uncorrected), with titles

    Parameters
    ----------
    version_labels : dict
        Mapping from leak-corrected version strings to human-readable labels.
        E.g., {"SP_v1.4.5_leak_corr": "Initial", ...}
    fiducial_version : str
        Fiducial version string (e.g., "SP_v1.4.11.2_leak_corr").

    Yields
    ------
    dict
        Figure specification with keys:
        - filename: output filename (e.g., "figure_v1.4.5.png")
        - version_leak_corr: leak-corrected version string
        - version_uncorr: uncorrected version string
        - leak_corrected: bool
        - title: title string or None (for paper figure)
        - is_paper_figure: bool
    """
    # Paper figure: fiducial, leak-corrected, no title
    yield {
        "filename": "figure.png",
        "version_leak_corr": fiducial_version,
        "version_uncorr": fiducial_version.replace("_leak_corr", ""),
        "leak_corrected": True,
        "title": None,
        "is_paper_figure": True,
    }

    # Dashboard figures: all versions × (corrected + uncorrected)
    # Sort by version string length descending to ensure proper matching
    for version_lc in sorted(version_labels.keys(), key=lambda v: -len(v)):
        version_num = extract_version_number(version_lc)
        label = version_labels[version_lc]
        version_uncorr = version_lc.replace("_leak_corr", "")

        # Leak-corrected figure
        yield {
            "filename": f"figure_{version_num}.png",
            "version_leak_corr": version_lc,
            "version_uncorr": version_uncorr,
            "leak_corrected": True,
            "title": label,
            "is_paper_figure": False,
        }

        # Uncorrected figure
        yield {
            "filename": f"figure_{version_num}_uncorrected.png",
            "version_leak_corr": version_lc,
            "version_uncorr": version_uncorr,
            "leak_corrected": False,
            "title": f"{label} (uncorrected)",
            "is_paper_figure": False,
        }


def get_version_figure_outputs(version_labels):
    """Get list of all per-version figure filenames for Snakemake outputs.

    Parameters
    ----------
    version_labels : dict
        Mapping from leak-corrected version strings to human-readable labels.

    Returns
    -------
    list of str
        List of 9 filenames (1 paper + 8 dashboard figures).
    """
    filenames = ["figure.png"]
    for version_lc in sorted(version_labels.keys(), key=lambda v: -len(v)):
        version_num = extract_version_number(version_lc)
        filenames.append(f"figure_{version_num}.png")
        filenames.append(f"figure_{version_num}_uncorrected.png")
    return filenames
