"""The pure-E/B decision variants, as one table every scan stage imports.

Three ASTRA decisions reach the pure-mode transform:
``pure_eb_integration_method`` (adaptive_dqags / fixed_operator),
``pure_eb_rebinning_weights`` (uniform / npairs / none) and
``pure_eb_integration_limits`` (fine_grid_extent / reporting_extent).  They do
not span a free cube, and the reason is worth stating precisely, because the
default table in ``pure_eb_transform`` confounds two effects that the scan
exists to separate.

The transform has two *geometries*:

* **pointwise** — evaluate the transform at the 20 reporting abscissae, take
  [tmin, tmax] from the reporting grid, apply no fine-to-reporting average
  (``rebinning_weights`` is ``none``).  The published chain.
* **averaged** — evaluate at the fine-grid nodes carrying reporting-bin
  support, take [tmin, tmax] from the fine-grid extent, and average the pure
  modes into the reporting bins under ``uniform`` or ``npairs`` weights.

``GEOMETRIES`` in ``pure_eb_transform`` maps each *method* to a default
geometry (dqags -> pointwise, fixed_operator -> averaged).  That is a table of
defaults, not a constraint, and taking it for a constraint is what makes the
published-vs-adopted comparison confounded: under the defaults, "not
pointwise" always also means "not dqags", so the paper's small-theta pure-B
excess cannot be attributed to the quadrature rather than the geometry.  The
two cross terms below -- dqags evaluated on the fine grid and averaged -- are
what break that confounding.

The reverse pairing -- ``fixed_operator`` in the pointwise geometry -- is
reachable only by approximation, and this table holds two readings of it.
``get_pure_EB_operator`` builds the single-grid operator, so evaluation points
must be nodes of the integration grid, and the 20 reporting ``meanr`` values
are up to 0.31% off the fine nodes; ``snap_eval=True`` evaluates at the nearest
node instead and records the offsets.  Snapping moves *both* terms of the
transform, and the two are not equally harmless:

* ``pointwise_fixedop`` takes the local ``0.5 * (xi_p +/- xi_m)`` term from the
  snapped fine node, i.e. from ~1/50 of a reporting bin's pairs.  That, not the
  quadrature, dominates its per-point noise.
* ``pointwise_fixedop_reporting`` snaps identically for the *integral* term but
  reads the local term from the 20-bin reporting measurement, exactly as the
  published dqags chain does (``local_from_eval_grid=True``).  It is the
  like-for-like operator counterpart of ``pointwise_dqags``.

That local-vs-integral split is a third axis, and it is not confined to the
pointwise geometry.  The averaged geometry also builds its local term from the
fine grid -- the block's own weighted mean of ``0.5 * (xi_p +/- xi_m)`` over the
reporting bin's fine nodes -- which is not the same number as the 20-bin
reporting measurement, because the latter is its own pair-weighted TreeCorr
measurement rather than a rebinning.  ``averaged_uniform_fixedop_reporting``
substitutes it, completing the square on that axis:

    local term from ...    pointwise geometry            averaged geometry
    the fine grid          pointwise_fixedop             averaged_uniform_fixedop
    the reporting xi       pointwise_fixedop_reporting   averaged_uniform_fixedop_reporting

``DROPPED`` holds the one genuine omission, ``averaged_npairs_dqags``, which is
merely expensive: seven variants run.

Option values are literals here rather than ``{decisions.*}`` in the recipes.
That is the point of a scan: its content must be the same in every universe.
A recipe reading the decisions would collapse to whichever single option the
active universe selects, and the comparison would vanish.  The decisions stay
decisions in the spec, where the ruling lives; this table measures them.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Variant:
    """One coherent (method, geometry, weights, limits) point of the scan."""

    name: str
    integration_method: str
    geometry: str
    rebinning_weights: str
    integration_limits: str
    label: str
    note: str
    #: Source the operator's local (non-integral) term from the *reporting*
    #: xi at the requested abscissa rather than from the fine grid at the
    #: snapped node.  Only meaningful for the snapped pointwise operator
    #: chain: dqags already reads the local term from the evaluation-grid xi,
    #: and the averaged geometry snaps nothing.  Default False, so every
    #: pre-existing variant is bit-identical.
    local_from_eval_grid: bool = False
    #: Name of the variant whose covariance bundle this one reads instead of
    #: building its own.  ``None`` (the default) means "own bundle".  Set only
    #: where the two propagated covariances are provably the same object, so
    #: that remaking one would be a fake dependency.  NOTE the borrowed bundle
    #: also carries the DONOR's data pure-mode vectors, so an aliased variant
    #: has no data vector of its own in the bundle and therefore no data PTE;
    #: pure_eb_scan_statistics reports that as absent rather than substituting
    #: the donor's number.
    covariance_from: str | None = None

    @property
    def covariance_source(self) -> str:
        """The variant name whose covariance bundle this one uses."""

        return self.covariance_from or self.name

    @property
    def transform_averaged(self) -> bool:
        return self.geometry == "averaged"

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "local_from_eval_grid": self.local_from_eval_grid,
            "covariance_from": self.covariance_from,
            "integration_method": self.integration_method,
            "geometry": self.geometry,
            "rebinning_weights": self.rebinning_weights,
            "integration_limits": self.integration_limits,
            "label": self.label,
            "note": self.note,
        }


#: The published chain, frozen in the ``archive_v2`` universe.
PUBLISHED = "pointwise_dqags"
#: The adopted chain (ruled 2026-08-31), the ``baseline`` universe.
ADOPTED = "averaged_uniform_fixedop"

VARIANTS: tuple[Variant, ...] = (
    Variant(
        "pointwise_dqags",
        "adaptive_dqags",
        "pointwise",
        "none",
        "reporting_extent",
        r"dqags, pointwise",
        "The published paper (archive_v2): adaptive QUADPACK at the 20 "
        "reporting abscissae, limits from the reporting grid.",
    ),
    Variant(
        "averaged_uniform_dqags",
        "adaptive_dqags",
        "averaged",
        "uniform",
        "fine_grid_extent",
        r"dqags, averaged (uniform)",
        "The cross term, and the reason the scan exists. With "
        "averaged_uniform_fixedop it isolates the quadrature effect at fixed "
        "geometry; with pointwise_dqags it isolates the geometry effect at "
        "fixed quadrature. Those two contrasts are what the default "
        "method->geometry mapping confounds.",
    ),
    Variant(
        "pointwise_fixedop",
        "fixed_operator",
        "pointwise",
        "none",
        "reporting_extent",
        r"operator, pointwise (nearest node)",
        "The fourth corner, reachable only by approximation. The single-grid "
        "operator evaluates at nodes, so each of the 20 reporting abscissae "
        "is snapped to the nearest fine node in log theta -- a shift of up to "
        "0.32%, recorded per point in snap_offsets and carried into the mode "
        "provenance. Completes the method x geometry square, so the "
        "quadrature and geometry effects can each be read along two edges "
        "rather than one. Read it alongside its offsets: the snap moves theta "
        "by the same order as the small-theta structure being measured.",
    ),
    Variant(
        "pointwise_fixedop_reporting",
        "fixed_operator",
        "pointwise",
        "none",
        "reporting_extent",
        r"operator, pointwise (reporting local term)",
        "The hybrid, and the corner pointwise_fixedop should have been. The "
        "pure-mode transform is a local term plus an integral term; "
        "pointwise_dqags takes the local term from the 20-bin reporting "
        "measurement and only the integral from the fine grid, while "
        "pointwise_fixedop -- forced to evaluate at a node -- takes BOTH from "
        "the snapped fine node, so its local term carries the shot noise of "
        "~1/50 of a reporting bin's pairs. This variant snaps exactly as "
        "pointwise_fixedop does for the integral term (same nearest-node "
        "choice, same recorded offsets) but restores the reporting-grid local "
        "term. It is therefore the like-for-like operator counterpart of the "
        "published pointwise chain: the difference from pointwise_dqags is "
        "quadrature (plus the snapped integral's theta shift), not fine-bin "
        "noise.",
        local_from_eval_grid=True,
    ),
    Variant(
        "averaged_uniform_fixedop",
        "fixed_operator",
        "averaged",
        "uniform",
        "fine_grid_extent",
        r"operator, averaged (uniform)",
        "ADOPTED (baseline): the fixed-quadrature linear operator, pure modes "
        "averaged into the reporting bins with equal weight.",
    ),
    Variant(
        "averaged_uniform_fixedop_reporting",
        "fixed_operator",
        "averaged",
        "uniform",
        "fine_grid_extent",
        r"operator, averaged (uniform), reporting local term",
        "The averaged geometry's counterpart of pointwise_fixedop_reporting, "
        "and the fourth reading of the local-term axis. Integral term exactly "
        "averaged_uniform_fixedop's -- the operator at every fine node "
        "carrying the bin's support, uniformly averaged -- but the local "
        "0.5*(xi_p +/- xi_m) term is the 20-bin reporting measurement rather "
        "than the uniform mean of the fine xi over the bin. Those two differ "
        "because the reporting xi is its own TreeCorr measurement (pair "
        "weighted within the bin), not a rebinning of the fine one. Note the "
        "MC covariance chains construct their reporting draws BY uniformly "
        "rebinning the fine draws, so there the two local terms coincide by "
        "construction and this variant's covariance matches "
        "averaged_uniform_fixedop's; the difference lives in the data and "
        "mock DATA VECTORS. That is why covariance_from aliases it rather "
        "than rebuilding: propagating 500 draws through a transform that "
        "provably returns the same covariance is a fake dependency (Cail, "
        "2026-09-06). Measured: substituting the reporting local term on a "
        "uniformly-rebinned draw moves the modes by 1.4e-20, against 1.8e-06 "
        "on the real reporting measurement. The cost is that the borrowed "
        "bundle carries the DONOR's data vector, so this variant has no data "
        "PTE; its signal is the mock <xi^B> shift, which comes from "
        "mock_pure_eb_modes_scan and is unaffected.",
        local_from_eval_grid=True,
        covariance_from="averaged_uniform_fixedop",
    ),
    Variant(
        "averaged_npairs_fixedop",
        "fixed_operator",
        "averaged",
        "npairs",
        "fine_grid_extent",
        r"operator, averaged (npairs)",
        "The ruled arm (fixed_operator universe): the same operator and "
        "window, averaged with TreeCorr pair counts (Asgari et al. 2019, "
        "Appendix A).",
    ),
)

#: Combinations deliberately absent, and why.  Kept in code so the scan's
#: coverage is self-documenting rather than an omission a reader must notice.
DROPPED = {
    "averaged_npairs_dqags": (
        "cost; the weights question is answered under the operator. Running "
        "the adaptive quadrature at ~1000 fine-grid nodes costs ~4x the "
        "operator per draw, and this variant would have bought only a second "
        "reading of the uniform-vs-npairs contrast that "
        "averaged_uniform_fixedop and averaged_npairs_fixedop already give at "
        "a quarter the price. Dropped by Cail 2026-09-04. Nothing structural "
        "prevents it: add it back to VARIANTS and it runs."
    ),
}

BY_NAME = {variant.name: variant for variant in VARIANTS}
NAMES = tuple(variant.name for variant in VARIANTS)


def resolve(names=None) -> list[Variant]:
    """Return the named variants in table order, or all of them."""

    if not names:
        return list(VARIANTS)
    unknown = [name for name in names if name not in BY_NAME]
    if unknown:
        dropped = [name for name in unknown if name in DROPPED]
        detail = ""
        if dropped:
            detail = " " + " ".join(f"{name}: {DROPPED[name]}" for name in dropped)
        raise ValueError(
            f"unknown variant(s) {unknown}; expected some of {list(NAMES)}.{detail}"
        )
    # Table order, not command-line order, so every stage's output is keyed
    # the same way regardless of how the recipe happened to list them.
    return [variant for variant in VARIANTS if variant.name in set(names)]


def add_variants_argument(parser, *, help_suffix=""):
    """Add the shared ``--variants`` flag."""

    parser.add_argument(
        "--variants",
        nargs="+",
        default=list(NAMES),
        choices=list(NAMES),
        metavar="NAME",
        help=(
            "Which decision variants to run, from the scan table "
            f"({', '.join(NAMES)}). Defaults to all seven." + help_suffix
        ),
    )
