"""Save a paper figure so its PDF is as wide as the paper prints it.

``bbox_inches="tight"`` crops (or grows) the canvas to the drawn artists, so a
figure created at the print width saves narrower or wider than that, and
LaTeX's ``width=`` rescales it -- fonts included.  ``savefig_at_width`` scales
the canvas (both dimensions, so the aspect is kept) until the cropped width is
the printed one; the style's font sizes are then the sizes on the page.
"""

import matplotlib.pyplot as plt


def savefig_at_width(fig, path, width, tol=1e-3, max_iter=10, **kwargs):
    """Save ``fig`` with a tight bounding box exactly ``width`` inches wide."""
    kwargs.setdefault("bbox_inches", "tight")
    pad = kwargs.get("pad_inches", plt.rcParams["savefig.pad_inches"])
    if pad == "layout":
        pad = 0.0
    for _ in range(max_iter):
        fig.canvas.draw()
        saved = fig.get_tightbbox(fig.canvas.get_renderer()).width + 2 * pad
        if abs(saved / width - 1.0) < tol:
            break
        fig.set_size_inches(fig.get_size_inches() * (width / saved))
    fig.savefig(path, **kwargs)
