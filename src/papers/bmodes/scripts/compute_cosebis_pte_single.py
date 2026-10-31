# %%
"""Compute the COSEBIS B-mode PTE matrix over the full scale-cut pair grid.

Reconciliation of the ASTRA ``cosebis_pte_per_cut`` output: the Snakemake DAG
*scattered* one JSON per (version, blind, i_min, i_max) cut via this script and
gathered them downstream. The spec describes the gathered NPZ that scatter never
materialised. Here a single CLI loops the same ``_pte_scale_cut_pairs()`` grid,
runs the *identical* per-pair COSEBI computation (same theta grid, nmodes, NaN-on-
convergence-failure handling), and writes one gathered NPZ — no scatter, no
per-pair deps to wire.
"""

import argparse
import os
import time

# Pin every BLAS/OpenMP/numba pool to one thread BEFORE numpy (and hence
# OpenBLAS) is imported. The pair loop below is parallelised across *processes*;
# a nested per-process thread pool would multiply into the node's 1280-process
# RLIMIT_NPROC and oversubscribe the SLURM cgroup (these libraries size their
# pools from hardware, not cgroup, concurrency). Module top rather than a pool
# initializer: the pool uses the *fork* context, so a child inherits the
# parent's already-sized pools and an initializer would run far too late. Fork
# is safe here because the parent never executes a numba parallel region before
# forking (numba's threading layer is not fork-safe once its threads exist) --
# calculate_cosebis is only ever called in a child. Assigned, not setdefault:
# an inherited value would silently defeat the pinning.
# The runner communicates the task's width via OMP_NUM_THREADS; capture it
# BEFORE the pinning below overwrites it, or the process pool (which sizes
# itself from this value at line ~206) silently collapses to width 1.
_POOL_WIDTH = int(os.environ.get("OMP_NUM_THREADS", "1"))

for _v in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "DUCC0_NUM_THREADS",
    "NUMBA_NUM_THREADS",
):
    os.environ[_v] = "1"

import multiprocessing  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402

import numpy as np  # noqa: E402
import treecorr  # noqa: E402
import tqdm  # noqa: E402
from plotting_utils import compute_chi2_pte  # noqa: E402

from sp_validation.b_modes import calculate_cosebis  # noqa: E402
from paper_config import load_config


def _pte_scale_cut_pairs(row=None):
    """(i_min, i_max) index pairs for the PTE matrix, excluding the polynomial-
    root-unstable subsets. Mirrors _pte_scale_cut_pairs() in claims.smk.

    ``row`` restricts the list to a single i_min row of the matrix. The full list
    runs i outermost, so a row is a contiguous subsequence of it: a row member's
    pairs keep their global relative order, and merge_cosebis_pte_rows.py can
    place each one back by its (i_min, i_max) key without re-sorting.
    """
    unstable = {(9, 10), (10, 11), (11, 12), (13, 14)}
    return [
        (i, j)
        for i in range(20)
        for j in range(i + 1, 21)
        if (i, j) not in unstable and (row is None or i == row)
    ]


def _compute_pair(gg, cov_path, nmodes, theta_min, theta_max):
    """COSEBI PTEs for one (theta_min, theta_max) cut.

    Bit-identical to the original per-pair scatter job. Returns a result dict,
    or None on a polynomial-root (sympy NoConvergence) failure.
    """
    try:
        results = calculate_cosebis(
            gg,
            nmodes=nmodes,
            scale_cuts=[(theta_min, theta_max)],
            cov_path=cov_path,
        )
    except Exception as e:
        if "NoConvergence" in str(type(e).__name__) or "NoConvergence" in str(e):
            return None
        raise

    result = results[(theta_min, theta_max)]
    En_full = result["En"]
    Bn_full = result["Bn"]
    cov_full = result["cov"]

    # 6 modes: first 6
    cov_B_6 = cov_full[nmodes : nmodes + 6, nmodes : nmodes + 6]
    cov_E_6 = cov_full[:6, :6]
    chi2_B_6, pte_B_6, _ = compute_chi2_pte(Bn_full[:6], cov_B_6)
    chi2_E_6, pte_E_6, _ = compute_chi2_pte(En_full[:6], cov_E_6)

    # 20 modes: all
    cov_B_20 = cov_full[nmodes:, nmodes:]
    cov_E_20 = cov_full[:nmodes, :nmodes]
    chi2_B_20, pte_B_20, _ = compute_chi2_pte(Bn_full, cov_B_20)
    chi2_E_20, pte_E_20, _ = compute_chi2_pte(En_full, cov_E_20)

    return {
        "En": En_full,
        "Bn": Bn_full,
        "chi2_B_6": chi2_B_6,
        "pte_B_6": pte_B_6,
        "chi2_E_6": chi2_E_6,
        "pte_E_6": pte_E_6,
        "chi2_B_20": chi2_B_20,
        "pte_B_20": pte_B_20,
        "chi2_E_20": chi2_E_20,
        "pte_E_20": pte_E_20,
    }


# Per-worker read-only context. Populated in the parent before the pool forks,
# so children inherit it (gg, the cov path and nmodes are identical for every
# pair) instead of re-pickling it per task.
_CTX = {}


def _worker(task):
    """Run one pair in a child process. Returns (k, result-or-None)."""
    k, theta_min, theta_max = task
    res = _compute_pair(
        _CTX["gg"], _CTX["cov_path"], _CTX["nmodes"], theta_min, theta_max
    )
    return k, res


def main(
    config,
    xi_integration,
    cov_integration,
    out_file,
    version=None,
    blind=None,
    row=None,
    jobs=None,
):
    t_start = time.time()
    fid = config["fiducial"]
    # Version/blind tag only the output filename + provenance record; the theta
    # grid, nmodes and per-pair COSEBI compute are version-independent (they read
    # from config["fiducial"]), so a sweep call over a non-fiducial catalog stays
    # bit-identical to the fiducial call save for the xi/cov inputs and the tag.
    version = version if version is not None else fid["version"]
    blind = blind if blind is not None else fid["blind"]
    nmodes = int(fid["nmodes"])  # 20 for full computation

    min_sep_int = fid["min_sep_int"]
    max_sep_int = fid["max_sep_int"]
    nbins_int = fid["nbins_int"]

    gg = treecorr.GGCorrelation(
        min_sep=min_sep_int, max_sep=max_sep_int, nbins=nbins_int, sep_units="arcmin"
    )
    gg.read(xi_integration)

    # Reporting theta grid (nbins+1 = 21 edges): geomspace(1', 250', 21)
    theta_grid = np.geomspace(fid["min_sep"], fid["max_sep"], fid["nbins"] + 1)
    n_grid = fid["nbins"] + 1

    # row is None for the whole grid; the cosebis_pte_grid sweep passes one
    # i_min per member and merge_cosebis_pte_rows.py reassembles the full NPZ.
    pairs = _pte_scale_cut_pairs(row=row)
    n_pairs = len(pairs)

    def _mat():
        return np.full((n_grid, n_grid), np.nan)

    pte_B_6_mat, chi2_B_6_mat = _mat(), _mat()
    pte_E_6_mat, chi2_E_6_mat = _mat(), _mat()
    pte_B_20_mat, chi2_B_20_mat = _mat(), _mat()
    pte_E_20_mat, chi2_E_20_mat = _mat(), _mat()

    pair_arr = np.array(pairs, dtype=int)
    theta_min_arr = np.full(n_pairs, np.nan)
    theta_max_arr = np.full(n_pairs, np.nan)
    En_arr = np.full((n_pairs, nmodes), np.nan)
    Bn_arr = np.full((n_pairs, nmodes), np.nan)
    convergence_failure = np.zeros(n_pairs, dtype=bool)

    tasks = []
    for k, (i_min, i_max) in enumerate(pairs):
        theta_min = theta_grid[i_min]
        theta_max = theta_grid[i_max]
        theta_min_arr[k] = theta_min
        theta_max_arr[k] = theta_max
        tasks.append((k, theta_min, theta_max))

    # Pairs are independent (_compute_pair reads only its arguments), so the
    # loop is a plain process pool. Results are stored by task index k, never by
    # completion order, so the output arrays are identical to the serial run.
    # Default to the runner's per-task width, NOT the machine's. sched_getaffinity
    # would report the whole node: lc runs every task as a subprocess of one dask
    # worker in one cgroup, and consumes no part of ASTRA's `resources:` block, so
    # nothing narrows a task's affinity to its slot. OMP_NUM_THREADS is the knob the
    # runner already sets for exactly this (see run_2pcf.py, run_pure_eb_semianalytic.sh).
    # Reading it from the environment rather than taking it as a recipe argument is
    # deliberate: how wide a task runs is a scheduling choice, not part of what it
    # computed, and lc keys staleness on the recipe command string -- so a core count
    # baked into the command would restale a scientifically identical artifact every
    # time the allocation changed shape.
    if jobs is None:
        jobs = _POOL_WIDTH
    n_jobs = jobs
    n_jobs = max(1, min(n_jobs, n_pairs))
    _CTX.update(gg=gg, cov_path=cov_integration, nmodes=nmodes)
    print(f"Running {n_pairs} scale-cut pairs over {n_jobs} processes")

    with ProcessPoolExecutor(
        max_workers=n_jobs, mp_context=multiprocessing.get_context("fork")
    ) as pool:
        futures = [pool.submit(_worker, t) for t in tasks]
        for fut in tqdm.tqdm(
            as_completed(futures), total=n_pairs, desc="COSEBIs scale cuts"
        ):
            k, res = fut.result()
            i_min, i_max = pairs[k]
            if res is None:
                convergence_failure[k] = True
                print(
                    f"WARNING: polyroots convergence failure for ({i_min}, {i_max}); NaN"
                )
                continue

            En_arr[k, :] = res["En"]
            Bn_arr[k, :] = res["Bn"]
            pte_B_6_mat[i_min, i_max] = res["pte_B_6"]
            chi2_B_6_mat[i_min, i_max] = res["chi2_B_6"]
            pte_E_6_mat[i_min, i_max] = res["pte_E_6"]
            chi2_E_6_mat[i_min, i_max] = res["chi2_E_6"]
            pte_B_20_mat[i_min, i_max] = res["pte_B_20"]
            chi2_B_20_mat[i_min, i_max] = res["chi2_B_20"]
            pte_E_20_mat[i_min, i_max] = res["pte_E_20"]
            chi2_E_20_mat[i_min, i_max] = res["chi2_E_20"]

    with open(out_file, "wb") as _f:
        np.savez(
            _f,
            version=version,
            blind=blind,
            nmodes=nmodes,
            mode_subsets=np.array([6, 20]),
            theta_grid=theta_grid,
            pairs=pair_arr,
            theta_min=theta_min_arr,
            theta_max=theta_max_arr,
            En=En_arr,
            Bn=Bn_arr,
            convergence_failure=convergence_failure,
            pte_B_6=pte_B_6_mat,
            chi2_B_6=chi2_B_6_mat,
            pte_E_6=pte_E_6_mat,
            chi2_E_6=chi2_E_6_mat,
            pte_B_20=pte_B_20_mat,
            chi2_B_20=chi2_B_20_mat,
            pte_E_20=pte_E_20_mat,
            chi2_E_20=chi2_E_20_mat,
        )
    print(f"Saved gathered COSEBI PTE matrices to {out_file}")
    print(f"[TIMING] Total ({n_pairs} pairs): {time.time() - t_start:.2f}s")


def _from_cli(argv=None):
    import yaml

    ap = argparse.ArgumentParser(
        description="Gathered COSEBI B-mode PTE matrix over the scale-cut pair grid (fiducial version + blind)."
    )
    ap.add_argument(
        "--config", required=True, help="Absolute path to bmodes config.yaml"
    )
    ap.add_argument(
        "--xi-integration",
        required=True,
        help="Fiducial 1000-bin integration-grid TreeCorr xi_pm .txt dump",
    )
    ap.add_argument(
        "--cov-integration",
        required=True,
        help="Fiducial 1000-bin Gaussian covariance (processed .txt)",
    )
    ap.add_argument(
        "--out-file",
        required=True,
        help="Output file (lc {output}): the single gathered PTE-matrix NPZ",
    )
    ap.add_argument(
        "--version",
        default=None,
        help="Catalog version tag (default: config.fiducial.version). Overridden "
        "by the version sweep to name the non-fiducial output NPZ.",
    )
    ap.add_argument(
        "--blind",
        default=None,
        help="Blind tag (default: config.fiducial.blind)",
    )
    ap.add_argument(
        "--row",
        type=int,
        default=None,
        help="Restrict the scale-cut pair grid to this i_min row (0-19), writing "
        "a partial NPZ of the same shape. Used by the cosebis_pte_grid row sweep; "
        "merge_cosebis_pte_rows.py gathers the rows back into the full-grid NPZ.",
    )
    ap.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Worker processes for the scale-cut pair loop. The default -- the "
        "number of CPUs in this process's affinity mask -- is the whole SLURM "
        "allocation, which is far too many when several members of the "
        "cosebis_pte_grid sweep run concurrently in one allocation. Callers that "
        "share an allocation MUST pass this explicitly; see that sweep's "
        "declaration for the arithmetic.",
    )
    a = ap.parse_args(argv)
    config = load_config(a.config)
    main(
        config,
        a.xi_integration,
        a.cov_integration,
        a.out_file,
        version=a.version,
        blind=a.blind,
        row=a.row,
        jobs=a.jobs,
    )


if __name__ == "__main__":
    _from_cli()
