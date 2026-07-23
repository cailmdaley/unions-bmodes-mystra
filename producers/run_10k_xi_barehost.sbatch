#!/bin/bash
# ============================================================================
# 10,000-bin convergence xi (two_point.xi_integration_convergence) — BARE-HOST.
#
# WHY bare-host (Cail's steer, 2026-07-05): TreeCorr needs no container, and the
# container path is an MPI dead end — its mpi4py is a MPICH wheel (libmpi.so.12)
# while the container's mpiexec is OpenMPI 4.1.4 and the historical launcher is
# OpenMPI 5.0.5; multi-node containerized MPI has no working launcher (no
# ssh/orted in the image, host<->host ssh publickey-denied, host/container ABI
# mismatch). Bare-host dissolves it: the host ~/.local mpi4py 4.1.1 links the
# SLURM-integrated /softs/openmpi/5.0.5-slurm-CentOS8 (libmpi.so.40 + libpmix.so.2),
# so `mpiexec --map-by ppr:N:node` does multi-node natively via SLURM PMIx. This
# is exactly the historical `xi_highres` snakemake rule's launcher (see
# workflow/rules/twopoint.smk) and Axel Guinot's proven run_corr.sh.
#
# run_2pcf_highres.py falls back to its standalone cat_config resolver when the
# full sp_validation stack is absent (bare-host), so NO PYTHONPATH is set — the
# resolver returns byte-identical shear config (path + leak_corr columns).
#
# Usage:
#   SCALED SMOKE (proves multi-node COMM_WORLD + output, ~min):
#     sbatch --nodes=2 --nodelist=n32,n33 -t 20:00 --export=ALL,SMOKE=1 run_10k_xi_barehost.sbatch
#   PRODUCTION (10k grid, ~1.5-2 days at 180 cores):
#     sbatch run_10k_xi_barehost.sbatch
# ============================================================================
#SBATCH --job-name=tenk_xi
#SBATCH -p pscompl
#SBATCH -N3
#SBATCH --nodelist=n32,n33,n35
#SBATCH --exclusive
#SBATCH -t 5-00:00:00
#SBATCH -x n09,n17,n36,n29,n31,n34
#SBATCH -o /n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/Reproductions/UNIONS/pure_eb/_scratch/tenk_xi/tenk_xi.%j.log

set -euo pipefail

WT=/automnt/n17data/cdaley/unions/code/sp_validation.worktrees/repro-paper-ii-astra
PROJ=/n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/Reproductions/UNIONS/pure_eb
SCRIPT=$WT/workflow/scripts/run_2pcf_highres.py
PY=/softs/intelpython/intelpython3-2025.1.0_196/bin/python3.12
OMPI=/softs/openmpi/5.0.5-slurm-CentOS8

# Host OpenMPI 5.0.5-slurm (ABI-matched to the host ~/.local mpi4py 4.1.1)
export PATH=$OMPI/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:${LD_LIBRARY_PATH:-}
# treecorr num_threads is read from SLURM_CPUS_PER_TASK by the script; pin it so
# ranks don't each grab a whole node's cores (mpiexec, not srun, launches ranks).
export SLURM_CPUS_PER_TASK=12
export OMP_NUM_THREADS=12

NODES=${SLURM_JOB_NUM_NODES:-3}
if [ "${SMOKE:-0}" = "1" ]; then
  NBINS=100; NPATCH=8; PPN=2
  OUT=$PROJ/_scratch/tenk_xi/smoke_2node
else
  NBINS=10000; NPATCH=50; PPN=5     # 3 nodes x 5 ranks x 12 threads = 180 cores
  # STAGING, not the lc output dir: `lc run xi_integration_convergence` clears
  # its {output} dir before running the (copy) recipe (snakemake removes declared
  # output first), which would delete a grid written straight there. The recipe
  # copies from this staging dir — which lc never touches — into {output}.
  OUT=$PROJ/_scratch/tenk_xi/production_out
fi
NTASKS=$(( NODES * PPN ))
mkdir -p "$OUT"

echo "[tenk_xi] $( [ "${SMOKE:-0}" = "1" ] && echo SMOKE || echo PRODUCTION ): nbins=$NBINS npatch=$NPATCH nodes=$NODES ppr=$PPN ntasks=$NTASKS -> $OUT"
echo "[tenk_xi] host mpiexec: $(which mpiexec)"; mpiexec --version | head -1

# NO PYTHONPATH -> run_2pcf_highres uses its standalone cat_config resolver.
mpiexec --map-by ppr:${PPN}:node -n ${NTASKS} \
  "$PY" "$SCRIPT" \
  --config     $WT/papers/bmodes/config/config.yaml \
  --cat-config $WT/cosmo_val/cat_config.yaml \
  --version    SP_v1.4.6.3_leak_corr \
  --nbins      ${NBINS} --npatch ${NPATCH} \
  --min-sep    0.5 --max-sep 300.0 \
  --out        "$OUT"

echo "[tenk_xi] mpiexec rc=$?  output -> $OUT"
ls -la "$OUT"/*npatch=1.txt 2>/dev/null || echo "[tenk_xi] WARNING: expected txt not found"
