"""The pure-mode covariance, and the χ² correction it implies.

The averaged fixed-operator pure-E/B estimator is exactly linear in the fine
ξ_± vector: ``ξ^{E/B/amb} = M · [ξ_+; ξ_−]``. So its covariance has two
constructions, selected by the ``pure_eb_covariance_method`` decision:

``exact_linear``
    ``M C_ξ Mᵀ``, with M built by pushing unit vectors through the pipeline's
    own :class:`~pure_eb_transform.AveragedPureEBTransform`. Noise-free, so no
    finite-sample (Hartlap) factor enters any χ² built on it.

``monte_carlo``
    The empirical covariance of N Gaussian draws of ξ_± pushed through the
    same transform — ``M Ĉ_ξ Mᵀ`` with Ĉ_ξ the sample covariance of the draws.
    Its inverse carries the Hartlap (N − p − 2)/(N − 1) factor.

Every covariance artifact records which construction made it
(``covariance_method``, plus ``n_samples`` for the Monte Carlo);
:func:`hartlap_samples` reads that back, so each χ² consumer applies the
correction the artifact it was handed calls for. Artifacts that carry no
``covariance_method`` are Monte Carlo with :data:`MC_SAMPLES_WITHOUT_RECORD`
draws.
"""

from __future__ import annotations

import numpy as np

COVARIANCE_METHODS = ("exact_linear", "monte_carlo")

#: Draw count of a pure-mode covariance artifact that records no construction:
#: every such artifact is a 20 × 100 Monte Carlo.
MC_SAMPLES_WITHOUT_RECORD = 2000


def add_covariance_method_argument(parser):
    """Add the required ``--covariance-method`` flag (the ASTRA decision)."""
    parser.add_argument(
        "--covariance-method",
        required=True,
        choices=COVARIANCE_METHODS,
        help="pure_eb_covariance_method: exact M C Mᵀ or the Monte Carlo.",
    )


def operator_matrix(transform, n_fine):
    """The (6·n_report × 2·n_fine) matrix of an averaged pure-mode transform.

    Column j is the transform applied to the j-th unit vector of the stacked
    fine ``[ξ_+; ξ_−]``. Refuses a transform that drops a non-finite fine
    node on any unit vector, since M would then not represent it.
    """
    columns = []
    for j in range(2 * n_fine):
        unit = np.zeros(2 * n_fine)
        unit[j] = 1.0
        columns.append(np.concatenate(transform(unit[:n_fine], unit[n_fine:])))
        if np.any(transform.last_dropped):
            raise ValueError(f"transform dropped non-finite output on unit vector {j}")
    return np.column_stack(columns)


def exact_covariance(transform, cov_xi):
    """``M C_ξ Mᵀ`` for an averaged pure-mode transform and a fine ξ_± covariance."""
    n_fine = cov_xi.shape[0] // 2
    M = operator_matrix(transform, n_fine)
    return M @ cov_xi @ M.T


def covariance_record(method, n_samples=None):
    """The keys a covariance artifact stores to name its construction."""
    if method == "exact_linear":
        return {"covariance_method": method}
    if method == "monte_carlo":
        return {"covariance_method": method, "n_samples": int(n_samples)}
    raise ValueError(f"unknown covariance method {method!r}")


def hartlap_samples(artifact):
    """The Hartlap sample count for χ² on this artifact's covariance.

    ``None`` when the covariance is exact (no correction); the Monte Carlo
    draw count otherwise. *artifact* is an open ``np.load`` result or any
    mapping with its keys.
    """
    keys = artifact.files if hasattr(artifact, "files") else artifact.keys()
    method = str(artifact["covariance_method"]) if "covariance_method" in keys else "monte_carlo"
    if method == "exact_linear":
        return None
    if "n_samples" in keys:
        return int(artifact["n_samples"])
    return MC_SAMPLES_WITHOUT_RECORD


def hartlap_factor(n_samples, dof):
    """(N − p − 2)/(N − 1), or exactly 1 when *n_samples* is ``None``."""
    return 1.0 if n_samples is None else (n_samples - dof - 2) / (n_samples - 1)
