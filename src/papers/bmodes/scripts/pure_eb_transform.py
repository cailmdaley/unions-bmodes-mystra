"""Pure E/B transform, resolved under the ``pure_eb_integration_method`` decision.

The Schneider (2022) decomposition that ``cosmo_numba`` implements *intends* a
linear map: the degree-``k`` interpolant it lays through the integrand samples
is linear in those samples.  ``get_pure_EB_modes`` then hands that interpolant
to QUADPACK's adaptive ``dqags``, which chooses its subdivision from the
sampled values — so the realised map is weakly non-linear, and on noisy input
its sub-tolerance quadrature error correlates with the noise and rectifies into
a spurious signal.  ``get_pure_EB_operator`` replaces the adaptive step by a
fixed Gauss-Legendre rule per grid cell, giving a data-independent matrix ``M``
with ``out = M @ [xi_p ; xi_m]``.

``pure_eb_integration_method`` in the ASTRA spec selects between them; this
module is the single place the choice is resolved.  Both paths take the same
geometry and the same per-sample inputs, and both return the reference
function's six-tuple in its own order
``(xip_E, xim_E, xip_B, xim_B, xip_amb, xim_amb)``.

Under ``adaptive_dqags`` the call is byte-for-byte the one the published paper
made.  Under ``fixed_operator`` the operator is built once per geometry (it
depends only on the grid, the bounds and the kernels) and applied as a matvec
per sample.

Geometry follows the method
---------------------------
``get_pure_EB_operator`` builds the operator for the *single-grid* form of the
transform: the evaluation points must be nodes of the integration grid
(``theta_eval`` a subset of ``theta``), because the operator's local
``0.5 * (xi_p(t) +/- xi_m(t))`` term is a column of that same grid.  The
paper's pointwise geometry — evaluate on the 20-bin reporting grid while
integrating on the 1000-bin grid — cannot be expressed that way.  So
``fixed_operator`` *defaults* every consuming chain to the single-grid
**transform-averaged** geometry — evaluate on the fine-grid nodes that fall
inside the reporting-bin span, with ``[tmin, tmax]`` at the fine-grid extent,
then average the pure modes into the reporting bins.  :func:`geometry_for` is
the single place that mapping is made, and :class:`AveragedPureEBTransform`
implements the averaged geometry.  ``adaptive_dqags`` keeps the published
pointwise geometry unchanged.

That mapping is a **default, not a constraint**, and the two directions are
not symmetric.  :class:`AveragedPureEBTransform` is method-agnostic by
construction, so the averaged geometry runs under ``adaptive_dqags`` just as
well -- the combination that separates the quadrature effect from the
geometry effect, which the defaults confound.  Callers reach it by passing
``geometry="averaged"`` to :func:`geometry_for` (or ``transform_averaged=True``
to the chain scripts, which already carry the flag).  The reverse pairing,
``fixed_operator`` in the pointwise geometry, needs the evaluation points
snapped onto the integration grid: opt in with ``snap_eval=True`` and read the
recorded offsets.  See :func:`geometry_for` for both.

The *weights* of that fine-to-reporting average are a separate axis:
``pure_eb_rebinning_weights`` in the ASTRA spec, reaching this module as
``--rebinning-weights`` and selecting ``npairs`` (the adopted Asgari et al.
2019 Appendix A pair-count weighting, applied to the whole pure mode at each
fine node so the collapsed local term is the reporting-bin xi itself),
``uniform`` (the unweighted mean, baseline 2026-08-31 to 2026-09-07), or ``none``
(the pointwise chain, which applies no fine-to-reporting rule at all).  The
two axes are separate but not independent, and
:func:`check_rebinning_weights` is where that is enforced: ``none`` is the
only coherent weighting in the pointwise geometry and never a coherent one in
the averaged geometry, so a mismatched universe selection fails loudly instead
of quietly running the wrong estimator.

:class:`TwoGridOperatorUnsupported` therefore no longer fires in normal
operation — it remains as a guard against a caller that builds a pointwise
:class:`PureEBTransform` under ``fixed_operator``, so that the evaluation points
are never silently snapped onto the fine grid (a ~0.3% shift in theta, right
where the small-theta B-mode behaviour under study lives).  See
``implementation-notes.md``.
"""

import numpy as np
from scipy import sparse

METHODS = ("adaptive_dqags", "fixed_operator")

# ``get_pure_EB_modes`` returns this order; the operator names its matrices the
# same way but orders them differently, so the mapping is written out.
REFERENCE_ORDER = ("xip_E", "xim_E", "xip_B", "xim_B", "xip_amb", "xim_amb")


GEOMETRIES = {"adaptive_dqags": "pointwise", "fixed_operator": "averaged"}

#: The *local* (non-integral) term of each pure-mode block, as the coefficients
#: of xi_+ and xi_- **at the evaluation point**.  Both implementations carry it
#: explicitly and identically -- ``schneider2022_nb`` writes
#: ``xip_E[i] = 0.5 * (xip[i] + xim[i] + int_tmp) - ...`` and
#: ``xip_B[i] = 0.5 * (xip[i] - xim[i] - int_tmp) - ...`` (and the same pair for
#: ``xim_E`` / ``xim_B``), while ``schneider2022_operator`` adds the same
#: numbers onto the operator row as ``out[ii] += 0.5`` / ``out[ii + n] += 0.5``
#: (E) or ``-= 0.5`` (B), ``ii`` being the evaluation point's own column.  The
#: two ambiguous blocks are pure integrals and have no local term.  Written out
#: here because :class:`PureEBTransform` needs to *re-source* that term when
#: ``local_from_eval_grid`` is set.
LOCAL_TERM_COEFFS = {
    "xip_E": (0.5, 0.5),
    "xim_E": (0.5, 0.5),
    "xip_B": (0.5, -0.5),
    "xim_B": (0.5, -0.5),
    "xip_amb": (0.0, 0.0),
    "xim_amb": (0.0, 0.0),
}

REBINNING_WEIGHTS = ("uniform", "npairs", "none")


def geometry_for(method, geometry=None):
    """Return the geometry a method implies, or an explicit *geometry* override.

    ``GEOMETRIES`` is a table of **defaults, not a constraint**.  The two axes
    are separable in one direction and genuinely bound in the other, and it is
    worth being precise about which is which:

    * **averaged + ``adaptive_dqags`` is well-defined and supported.**  Nothing
      in :class:`AveragedPureEBTransform` is operator-specific: it chooses the
      evaluation nodes, the window and the fine-to-reporting average, then
      delegates one transform call to :class:`PureEBTransform`, which
      dispatches on the method.  Handing it ``adaptive_dqags`` evaluates the
      published quadrature at the fine-grid abscissae and averages the result
      under either weighting.  This is the combination that separates the
      *quadrature* effect from the *geometry* effect, which the default table
      confounds — under the defaults alone, "not pointwise" always also means
      "not dqags".  It costs ~50x a matvec, because dqags runs per evaluation
      node and the fine grid has ~1000 of them rather than 20.
    * **pointwise + ``fixed_operator`` requires snapping, and says so.**
      ``get_pure_EB_operator`` builds the *single-grid* operator, whose local
      term is a column of the integration grid, so evaluation points must be
      nodes of that grid; the 20 reporting ``meanr`` values are all off the
      fine nodes.  Passing ``snap_eval=True`` to :class:`PureEBTransform`
      evaluates at the nearest node in *log* theta instead (the grid is
      geometric, so a linear metric would bias the choice outward), records
      the per-point offsets on ``snap_offsets`` / ``snap_offset_max``, and
      carries them into the mode provenance.  Without it the constructor
      still raises :class:`TwoGridOperatorUnsupported`, which is what keeps an
      averaged caller from silently getting a two-grid operator.

      The snap is a real approximation, not a formality: it moves theta by up
      to ~0.32%, the same order as the small-theta B-mode structure under
      study.  That is why it is opt-in and recorded rather than absorbed --
      the variant it enables is interpretable only alongside its offsets.

    Passing *geometry* explicitly is therefore how a caller says "I mean this
    geometry, whatever the method's default is".  Every existing call site
    passes one argument and is unaffected.
    """

    if method not in GEOMETRIES:
        raise ValueError(
            f"unknown integration method {method!r}; expected one of {METHODS}"
        )
    if geometry is None:
        return GEOMETRIES[method]
    if geometry not in ("pointwise", "averaged"):
        raise ValueError(
            f"unknown geometry {geometry!r}; expected 'pointwise' or 'averaged'"
        )
    return geometry


#: ``pure_eb_integration_limits`` options -> (theta_min, theta_max) transform
#: window in arcmin, applied by subsetting the fine integration grid before
#: the operator is built.  ``None`` bounds mean "the grid's own edge".
INTEGRATION_LIMITS = {
    # The data fine grid ends at 300'; the GLASS fine grid at 500'.  Capping
    # at 300' everywhere keeps mock and data operators on the same extent.
    # The wide floor is the default because the averaged operator needs
    # support below the 1' reporting floor at its edge.
    "fine_grid_extent": (None, 300.0),
    # The published pointwise chain takes [theta_min, theta_max] from the
    # reporting grid inside the transform itself, but its supplied fine grid
    # is capped at 300' to match the data-chain extent.  This is a no-op for
    # the data grid and preserves the historical cap for the 500' mock grid.
    # Pointwise-only.
    "reporting_extent": (None, 300.0),
}


def resolve_integration_limits(option, geometry):
    """Map a ``pure_eb_integration_limits`` option to a transform window.

    Returns ``(theta_min, theta_max)`` with ``None`` bounds meaning the
    supplied grid's own edge.  For pointwise geometry, the transform still
    takes ``tmin``/``tmax`` from the reporting grid internally, while the
    supplied fine arrays are capped at 300' to match the data-chain extent.
    Coherence with the geometry is enforced the same way
    :func:`check_rebinning_weights` does it: ``reporting_extent`` describes
    the pointwise chain only, and the averaged geometry must say which window
    it runs in.
    """

    if option not in INTEGRATION_LIMITS:
        raise ValueError(
            f"unknown integration limits {option!r}; expected one of "
            f"{sorted(INTEGRATION_LIMITS)}"
        )
    if geometry == "pointwise":
        if option != "reporting_extent":
            raise ValueError(
                "the pointwise geometry takes [tmin, tmax] from the reporting "
                "grid inside the transform; --integration-limits must be "
                f"'reporting_extent', not {option!r}."
            )
        return INTEGRATION_LIMITS[option]
    if option == "reporting_extent":
        raise ValueError(
            "'reporting_extent' describes the pointwise chain; the averaged "
            "geometry needs a fine-grid window: pass 'fine_grid_extent'."
        )
    return INTEGRATION_LIMITS[option]


def check_rebinning_weights(rebinning_weights, geometry):
    """Assert the fine-to-reporting weighting is coherent with the geometry.

    The averaged geometry always reduces 1000 fine nodes to 20 reporting
    bins, so it always needs a rule and ``none`` is never one.  The pointwise
    geometry evaluates the transform at the 20 reporting abscissae directly,
    so no fine-to-reporting rule is applied to the pure modes and ``none`` is
    the only coherent selection.  Checked rather than silently ignored: a
    universe that pairs the two wrongly should fail here, not run the other
    estimator.
    """

    if rebinning_weights not in REBINNING_WEIGHTS:
        raise ValueError(
            f"unknown rebinning weights {rebinning_weights!r}; expected one of "
            f"{REBINNING_WEIGHTS}"
        )
    if geometry == "averaged" and rebinning_weights == "none":
        raise ValueError(
            "the transform-averaged geometry averages the pure modes from the "
            "fine grid into the reporting bins, so it needs a rule: pass "
            "--rebinning-weights uniform or npairs. 'none' describes the "
            "pointwise chain, where no such average is taken."
        )
    if geometry == "pointwise" and rebinning_weights != "none":
        raise ValueError(
            "the pointwise geometry evaluates the transform at the 20 "
            "reporting abscissae and applies no fine-to-reporting rule to the "
            "pure modes, so --rebinning-weights must be 'none', not "
            f"{rebinning_weights!r}."
        )
    return rebinning_weights


class TwoGridOperatorUnsupported(NotImplementedError):
    """Raised when ``fixed_operator`` meets a two-grid (pointwise) call."""


def add_integration_method_argument(parser):
    """Add the required ``--integration-method`` flag.

    Deliberately has no default: the ASTRA spec is the single source of truth
    for which method a given universe adopts.
    """
    parser.add_argument(
        "--integration-method",
        required=True,
        choices=METHODS,
        help=(
            "Quadrature for the pure E/B integral transform: 'adaptive_dqags' "
            "is the published paper's adaptive QUADPACK path (non-linear on "
            "noisy input); 'fixed_operator' builds the fixed-quadrature linear "
            "operator M and applies it as a matvec."
        ),
    )


def add_integration_limits_argument(parser):
    """Add the required ``--integration-limits`` flag.

    Deliberately has no default, for the same reason the other decision flags
    have none: the ASTRA spec (``pure_eb_integration_limits``) is the single
    source of truth for the transform window a given universe adopts.
    """
    parser.add_argument(
        "--integration-limits",
        required=True,
        choices=sorted(INTEGRATION_LIMITS),
        help=(
            "[theta_min, theta_max] window of the pure-mode transform: "
            "'fine_grid_extent' keeps the grid up to 300' (default), and "
            "'reporting_extent' names the pointwise chain's own limits "
            "(pointwise-only, no fine-grid subsetting)."
        ),
    )


def add_rebinning_weights_argument(parser):
    """Add the required ``--rebinning-weights`` flag.

    Deliberately has no default, for the same reason
    :func:`add_integration_method_argument` has none: the ASTRA spec is the
    single source of truth for what a given universe adopts.
    """
    parser.add_argument(
        "--rebinning-weights",
        required=True,
        choices=REBINNING_WEIGHTS,
        help=(
            "Weights for the fine-to-reporting average of the pure modes: "
            "'npairs' is the adopted pair-count weighting of Asgari et al. "
            "(2019) Appendix A, 'uniform' the unweighted mean, and 'none' the "
            "pointwise chain, which takes no such average. Must match the "
            "geometry the integration method implies."
        ),
    )


class PureEBTransform:
    """One geometry's pure E/B transform, under the selected method.

    The geometry — evaluation grid, integration grid, ``[tmin, tmax]``, padding
    and interpolation order — is bound at construction; each call transforms
    one ``(xip, xim, xip_int, xim_int)`` sample.  Under ``fixed_operator`` the
    operator is built eagerly here, so an unsupported geometry fails before any
    Monte-Carlo draw is made.
    """

    def __init__(
        self,
        method,
        *,
        theta,
        theta_int,
        tmin,
        tmax,
        pad_xim=True,
        pad_theta_max_decade=1.0,
        interp_order=5,
        parallel=False,
        snap_eval=False,
        local_from_eval_grid=False,
    ):
        if method not in METHODS:
            raise ValueError(
                f"unknown integration method {method!r}; expected one of {METHODS}"
            )
        self.method = method
        self.theta = np.ascontiguousarray(theta, dtype=float)
        self.theta_int = np.ascontiguousarray(theta_int, dtype=float)
        self.tmin = float(tmin)
        self.tmax = float(tmax)
        self.pad_xim = bool(pad_xim)
        self.pad_theta_max_decade = float(pad_theta_max_decade)
        self.interp_order = int(interp_order)
        self.parallel = bool(parallel)
        self.snap_eval = bool(snap_eval)
        # Source the local (non-integral) term from the caller's evaluation-grid
        # xi rather than from the integration vector at the evaluation node.
        # A no-op unless the evaluation points were snapped -- without a snap
        # the two are the same numbers, and _apply_operator checks that they
        # are.  With a snap they are not: the local term would otherwise come
        # from a single fine node (~1/50 of a reporting bin's pairs), which is
        # what makes the snapped pointwise operator chain noisy.  Setting this
        # takes it from the reporting measurement instead, exactly as
        # get_pure_EB_modes does on the pointwise dqags path.
        self.local_from_eval_grid = bool(local_from_eval_grid)
        if self.local_from_eval_grid and method != "fixed_operator":
            raise ValueError(
                "local_from_eval_grid describes how the fixed_operator matvec "
                "sources its local term; adaptive_dqags already reads it from "
                "the evaluation-grid xi, so the flag is meaningless there."
            )
        # Per-evaluation-point |log(theta_snapped/theta_requested)| in dex, and
        # the same as a signed fractional offset. Empty unless snapping ran.
        self.snap_offsets = np.zeros(0, dtype=float)
        self.snap_offset_max = 0.0
        self._matrices = None
        if method == "fixed_operator":
            self._build_operator()

    # -- fixed_operator ----------------------------------------------------
    def _build_operator(self):
        from cosmo_numba.B_modes.schneider2022_operator import get_pure_EB_operator

        n_int = self.theta_int.size
        # Nearest node in LOG theta: the grid is geometric, so a linear metric
        # would bias the choice outward at every point.
        log_int = np.log(self.theta_int)
        log_eval = np.log(self.theta)
        idx = np.clip(np.searchsorted(log_int, log_eval), 1, n_int - 1)
        take_left = np.abs(log_int[idx - 1] - log_eval) <= np.abs(
            log_int[idx] - log_eval
        )
        nearest = np.where(take_left, idx - 1, idx)
        off_grid = int(
            np.count_nonzero(
                ~np.isclose(self.theta_int[nearest], self.theta, rtol=1e-12, atol=0.0)
            )
        )
        if off_grid and not self.snap_eval:
            raise TwoGridOperatorUnsupported(
                "fixed_operator needs the evaluation grid to be a subset of the "
                "integration grid: get_pure_EB_operator builds the single-grid "
                f"operator, and {off_grid} of {self.theta.size} evaluation "
                "points are not nodes of the integration grid. Pass "
                "snap_eval=True to evaluate at the nearest node instead (the "
                "pointwise fixed_operator chain), accepting the theta shift "
                "that implies; see implementation-notes.md."
            )
        self._eval_index = nearest
        if off_grid:
            # Snapping MOVES the evaluation points. On the fiducial pairing
            # (20 reporting meanr against the 1000-bin fine grid) the shift is
            # <= 0.32%, which is small but not nothing: it is the same order as
            # the small-theta B-mode structure this scan exists to measure, so
            # it is recorded rather than absorbed, and travels with the modes.
            self.snap_offsets = self.theta_int[nearest] / self.theta - 1.0
            self.snap_offset_max = float(np.max(np.abs(self.snap_offsets)))
            # Evaluate where the operator can actually be built.
            self.theta = np.ascontiguousarray(self.theta_int[nearest], dtype=float)

        operator = get_pure_EB_operator(
            theta_eval=self.theta,
            theta=self.theta_int,
            tmin=self.tmin,
            tmax=self.tmax,
            pad_xim=self.pad_xim,
            pad_theta_max_decade=self.pad_theta_max_decade,
            interp_order=self.interp_order,
            outputs=REFERENCE_ORDER,
        )
        # Rows whose t-dependent support holds fewer than interp_order + 1
        # nodes carry NaN and are flagged invalid. The fiducial evaluation
        # ranges sit strictly inside [tmin, tmax] with ample support, so an
        # invalid row means the geometry is not the one that was validated.
        invalid = {
            name: int(np.count_nonzero(~mask))
            for name, mask in operator["valid"].items()
            if not mask.all()
        }
        if invalid:
            raise ValueError(
                "fixed_operator produced under-determined (NaN) rows — the "
                f"evaluation grid is too close to [tmin, tmax]: {invalid}"
            )
        self._matrices = {
            name: np.ascontiguousarray(operator["matrices"][name])
            for name in REFERENCE_ORDER
        }

    def _apply_operator(self, xip, xim, xip_int, xim_int):
        # The operator's local term is a column of the integration grid, so the
        # evaluation-grid values must be the integration-grid values at the
        # shared nodes. That holds by construction for every single-grid
        # caller; a caller that broke it would get a silently different
        # estimator rather than an error, so it is checked.
        xip_int = np.asarray(xip_int, dtype=float)
        xim_int = np.asarray(xim_int, dtype=float)
        # When the evaluation points were snapped, the caller's own eval-grid
        # values are xi at the REQUESTED theta, which is deliberately not where
        # the operator evaluates. The local term must read the integration
        # vector at the snapped nodes, so the check is skipped rather than
        # failed -- the snap, not a mismatch, is what makes them differ.
        checks = (
            ()
            if self.snap_offsets.size
            else (
                ("xip", xip, xip_int[self._eval_index]),
                ("xim", xim, xim_int[self._eval_index]),
            )
        )
        for name, values, reference in checks:
            if not np.array_equal(np.asarray(values, dtype=float), reference):
                raise ValueError(
                    f"fixed_operator: {name} on the evaluation grid differs from "
                    f"{name}_int at the same nodes; the operator's local term "
                    "reads the integration vector, so the two must agree."
                )
        stacked = np.ascontiguousarray(np.concatenate([xip_int, xim_int]))
        modes = tuple(self._matrices[name] @ stacked for name in REFERENCE_ORDER)
        if not self.local_from_eval_grid:
            return modes
        # Re-source the local term.  The operator built it from the integration
        # vector at the evaluation node (LOCAL_TERM_COEFFS, added onto column
        # ``ii`` when the row was assembled), so subtracting that and adding the
        # same combination of the caller's evaluation-grid values leaves the
        # integral term untouched and swaps the local one.  Identically zero
        # when nothing was snapped, since then the two agree by construction.
        delta_p = np.asarray(xip, dtype=float) - xip_int[self._eval_index]
        delta_m = np.asarray(xim, dtype=float) - xim_int[self._eval_index]
        return tuple(
            values
            + LOCAL_TERM_COEFFS[name][0] * delta_p
            + LOCAL_TERM_COEFFS[name][1] * delta_m
            for name, values in zip(REFERENCE_ORDER, modes)
        )

    # -- adaptive_dqags ----------------------------------------------------
    def _apply_dqags(self, xip, xim, xip_int, xim_int):
        from cosmo_numba.B_modes.schneider2022 import get_pure_EB_modes

        return get_pure_EB_modes(
            theta=self.theta,
            xip=xip,
            xim=xim,
            theta_int=self.theta_int,
            xip_int=xip_int,
            xim_int=xim_int,
            tmin=self.tmin,
            tmax=self.tmax,
            parallel=self.parallel,
            pad_xim=self.pad_xim,
            pad_theta_max_decade=self.pad_theta_max_decade,
            interp_order=self.interp_order,
        )

    def __call__(self, xip, xim, xip_int, xim_int):
        """Return ``(xip_E, xim_E, xip_B, xim_B, xip_amb, xim_amb)``."""
        if self.method == "adaptive_dqags":
            return self._apply_dqags(xip, xim, xip_int, xim_int)
        return self._apply_operator(xip, xim, xip_int, xim_int)


def make_binning_matrix(
    left_edges, right_edges, theta_int, npairs=None, *, weights
):
    """Build the normalized fine-to-reporting matrix under *weights*.

    ``weights="npairs"`` gives each fine node the pair count TreeCorr
    measured at that separation; ``weights="uniform"`` gives every node in a
    reporting bin equal weight and reads ``npairs`` not at all, which is why
    the callers must not require the column under that selection.
    Normalisation is the same either way — each row sums to one.  So is the
    validity masking, up to the pair-count term it cannot have: ``npairs`` is
    strictly positive on the fiducial fine grid, so the two masks coincide
    there.
    """

    check_rebinning_weights(weights, "averaged")
    left_edges = np.asarray(left_edges, dtype=float)
    right_edges = np.asarray(right_edges, dtype=float)
    theta_int = np.asarray(theta_int, dtype=float)
    n_reporting = left_edges.size
    reporting_edges = np.concatenate([left_edges, [right_edges[-1]]])
    bin_indices = np.digitize(theta_int, reporting_edges) - 1
    valid = (bin_indices >= 0) & (bin_indices < n_reporting)
    if weights == "npairs":
        if npairs is None:
            raise ValueError("npairs weighting needs the integration npairs column")
        npairs = np.asarray(npairs, dtype=float)
        if npairs.shape != theta_int.shape:
            raise ValueError("npairs must have one entry per integration bin")
        if not np.all(np.isfinite(npairs)) or np.any(npairs < 0):
            raise ValueError("npairs must be finite and non-negative")
        valid &= npairs > 0
        weight_column = npairs[valid]
    else:
        weight_column = np.ones(int(np.count_nonzero(valid)), dtype=float)
    matrix = sparse.csr_matrix(
        (weight_column, (bin_indices[valid], np.where(valid)[0])),
        shape=(n_reporting, theta_int.size),
    )
    row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(row_sums == 0):
        raise ValueError("At least one reporting bin has no integration samples")
    return sparse.diags(1.0 / row_sums) @ matrix


def average_modes(modes, binning_matrix, fine_indices, companions=None):
    """Average finite fine-grid modes, renormalizing each block separately.

    *companions* is an optional pair of fine-grid vectors (in practice
    ``(xip, xim)`` on ``fine_indices``) averaged with **each block's own**
    normalized matrix.  The block matrices differ whenever a block dropped
    non-finite fine bins, and the caller that needs this -- the reporting
    local-term substitution -- must subtract exactly the average the block
    actually took, not a nominal one.  Returns ``(averaged, dropped)`` when
    *companions* is None and ``(averaged, dropped, per_block)`` otherwise,
    ``per_block`` mapping each block name to its averaged companions.
    """

    matrix = binning_matrix[:, fine_indices]
    averaged = []
    dropped = []
    per_block = {}
    for block_name, values in zip(REFERENCE_ORDER, modes):
        values = np.asarray(values, dtype=float)
        finite = np.isfinite(values)
        dropped.append(int(np.count_nonzero(~finite)))
        block_matrix = matrix[:, finite]
        row_sums = np.asarray(block_matrix.sum(axis=1)).ravel()
        if np.any(row_sums == 0):
            empty = np.flatnonzero(row_sums == 0).tolist()
            raise ValueError(
                f"All finite fine bins were dropped for {block_name} in "
                f"reporting bins {empty}"
            )
        block_matrix = sparse.diags(1.0 / row_sums) @ block_matrix
        averaged.append(np.asarray(block_matrix @ values[finite]).ravel())
        if companions is not None:
            per_block[block_name] = tuple(
                np.asarray(block_matrix @ np.asarray(c, dtype=float)[finite]).ravel()
                for c in companions
            )
    if companions is None:
        return tuple(averaged), np.asarray(dropped, dtype=int)
    return tuple(averaged), np.asarray(dropped, dtype=int), per_block


class AveragedPureEBTransform:
    """The transform-averaged (single-grid) geometry, under either method.

    The transform is evaluated on the fine-grid nodes that carry reporting-bin
    support and the resulting pure modes are averaged into the reporting bins
    under the ``rebinning_weights`` selection (``uniform`` or ``npairs``; see
    :func:`make_binning_matrix`).  ``[tmin, tmax]`` is the *fine-grid* extent
    rather than the reporting edges: evaluating at fine bins abutting the reporting
    limits leaves fewer than ``interp_order + 1`` quadrature nodes there
    (degenerate spline -> silent garbage), whereas the fine-grid limits keep
    real support at every evaluation point.  Note this defines the pure modes
    over the wider separation interval.

    Each call takes the fine-grid ``(xip_int, xim_int)`` and returns the six
    reporting-bin arrays in :data:`REFERENCE_ORDER`; the per-block count of
    dropped non-finite fine bins is left on :attr:`last_dropped`.

    ``theta_window=(lo, hi)`` (``None`` bounds = the grid's own edge) subsets
    the fine grid before the operator is built — the
    ``pure_eb_integration_limits`` decision, resolved by
    :func:`resolve_integration_limits`.  Calls still take *full-length*
    fine-grid inputs and the subsetting is applied internally, so every call
    site passes its arrays unchanged.
    """

    def __init__(
        self,
        method,
        *,
        left_edges,
        right_edges,
        theta_int,
        rebinning_weights,
        npairs=None,
        theta_window=None,
        local_from_eval_grid=False,
        **transform_kwargs,
    ):
        self.method = method
        self.rebinning_weights = rebinning_weights
        # Take the local 0.5*(xi_p +/- xi_m) term from the reporting-grid xi
        # rather than from this geometry's own fine-to-reporting average of it.
        # See __call__: the substitution happens AFTER the average, per block,
        # because each block's normalized matrix is what actually weighted the
        # local term into the reporting bin.
        self.local_from_eval_grid = bool(local_from_eval_grid)
        theta_int = np.ascontiguousarray(theta_int, dtype=float)
        self._n_input = len(theta_int)
        if theta_window is not None:
            lo, hi = theta_window
            keep = np.ones(self._n_input, dtype=bool)
            if lo is not None:
                keep &= theta_int >= float(lo)
            if hi is not None:
                keep &= theta_int <= float(hi)
            if np.count_nonzero(keep) < 10:
                raise ValueError(
                    f"theta_window {theta_window} leaves too few fine bins "
                    f"({np.count_nonzero(keep)})"
                )
            self._keep = keep
            theta_int = theta_int[keep]
            if npairs is not None:
                npairs = np.asarray(npairs, dtype=float)[keep]
        else:
            self._keep = None
        self.theta_int = np.ascontiguousarray(theta_int, dtype=float)
        self.binning_matrix = make_binning_matrix(
            left_edges,
            right_edges,
            self.theta_int,
            npairs,
            weights=rebinning_weights,
        )
        self.fine_indices = np.flatnonzero(
            np.asarray(self.binning_matrix.sum(axis=0)).ravel() > 0
        )
        self.tmin = float(self.theta_int.min()) * (1 - 1e-9)
        self.tmax = float(self.theta_int.max()) * (1 + 1e-9)
        self.transform = PureEBTransform(
            method,
            theta=self.theta_int[self.fine_indices],
            theta_int=self.theta_int,
            tmin=self.tmin,
            tmax=self.tmax,
            **transform_kwargs,
        )
        self.last_dropped = np.zeros(len(REFERENCE_ORDER), dtype=int)

    def __call__(self, xip_int, xim_int, xip=None, xim=None):
        """Return the six reporting-bin arrays in :data:`REFERENCE_ORDER`.

        *xip* / *xim* are the reporting-grid xi and are used only under
        ``local_from_eval_grid``, where they replace this geometry's own
        fine-to-reporting average of the local term.  Required in that case
        and ignored otherwise, which keeps every existing two-argument call
        site unchanged.
        """

        if self.local_from_eval_grid and (xip is None or xim is None):
            raise ValueError(
                "local_from_eval_grid needs the reporting-grid xi: call "
                "transform(xip_int, xim_int, xip, xim)."
            )
        xip_int = np.asarray(xip_int, dtype=float)
        xim_int = np.asarray(xim_int, dtype=float)
        if self._keep is not None:
            if len(xip_int) != self._n_input or len(xim_int) != self._n_input:
                raise ValueError(
                    f"expected full-grid inputs of length {self._n_input}, "
                    f"got {len(xip_int)}/{len(xim_int)}"
                )
            xip_int = xip_int[self._keep]
            xim_int = xim_int[self._keep]
        modes = self.transform(
            xip_int[self.fine_indices], xim_int[self.fine_indices], xip_int, xim_int
        )
        if not self.local_from_eval_grid:
            averaged, self.last_dropped = average_modes(
                modes, self.binning_matrix, self.fine_indices
            )
            return averaged
        # Each fine node's mode carries the local term at that node, so after
        # the average the reporting bin holds the block's own weighted mean of
        # 0.5*(xi_p +/- xi_m) over the bin's fine nodes.  Subtract exactly that
        # mean and add the reporting measurement in its place; the integral
        # term is untouched, and the delta is zero where the two coincide.
        averaged, self.last_dropped, local_means = average_modes(
            modes,
            self.binning_matrix,
            self.fine_indices,
            companions=(xip_int[self.fine_indices], xim_int[self.fine_indices]),
        )
        xip = np.asarray(xip, dtype=float)
        xim = np.asarray(xim, dtype=float)
        return tuple(
            values
            + LOCAL_TERM_COEFFS[name][0] * (xip - local_means[name][0])
            + LOCAL_TERM_COEFFS[name][1] * (xim - local_means[name][1])
            for name, values in zip(REFERENCE_ORDER, averaged)
        )
