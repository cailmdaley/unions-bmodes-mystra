"""Realisation-level process pool for the GLASS mock scatter scripts.

Both scatter transforms are a serial loop over 350 independent realisations
whose inner numba kernels parallelise poorly.  Spreading the realisations
over processes and pinning each worker to a single thread uses the
allocation far better; the pool is ordered (`Pool.map`), so the stacked
output is bit-identical to the serial loop.
"""

from __future__ import annotations

import multiprocessing
import os

_STATE: dict = {}


def pool_size(n_items: int) -> int:
    """Workers to use: the thread pin the runner sets, capped by the work available."""

    try:
        requested = int(os.environ.get("OMP_NUM_THREADS", "1"))
    except ValueError:
        requested = 1
    return max(1, min(requested, n_items))


def _initializer(state: dict) -> None:
    """Pin the worker to one thread and hand it the shared read-only state."""

    for var in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS"):
        os.environ[var] = "1"
    try:
        import numba

        # Adopt the 1-thread environment in numba's config *before* anything
        # launches the thread pool: set_num_threads() launches it, and a later
        # reload_config() (every JIT compile) raises if the cached
        # NUMBA_NUM_THREADS global still holds the parent's wider pin.
        try:
            numba.core.config.reload_config()
        except RuntimeError:  # threads already launched pre-fork; pin below
            pass
        numba.set_num_threads(1)
    except Exception:  # numba may be absent or already fixed at one thread
        pass
    _STATE.clear()
    _STATE.update(state)


def state() -> dict:
    """The state the parent handed the workers (or the parent's own, serially)."""

    return _STATE


def map_realisations(worker, indices, state_dict):
    """Run *worker* over *indices* in order, with *state_dict* shared by fork."""

    workers = pool_size(len(indices))
    if workers == 1:
        _initializer(state_dict)
        return [worker(index) for index in indices]
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(
        processes=workers, initializer=_initializer, initargs=(state_dict,)
    ) as pool:
        return pool.map(worker, indices, chunksize=1)
