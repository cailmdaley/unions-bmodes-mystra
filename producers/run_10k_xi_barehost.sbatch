#!/bin/bash
# The 10,000-bin xi_pm run behind data_vectors.xi_integration_convergence: TreeCorr
# on the leakage-corrected size-cut catalog over 0.5'-300', split across MPI
# ranks with npatch = 50. It runs outside the analysis's container because the
# measurement needs multi-node MPI, and it writes to the staging directory that
# the input `tenk_xi_staging` points at.
#
#   sbatch run_10k_xi_barehost.sbatch
#   sbatch --nodes=2 -t 20:00 --export=ALL,SMOKE=1 run_10k_xi_barehost.sbatch   # small test
#SBATCH --job-name=tenk_xi
#SBATCH -p pscompl
#SBATCH -N3
#SBATCH --nodelist=n32,n33,n35
#SBATCH --exclusive
#SBATCH -t 5-00:00:00
#SBATCH -x n09,n17,n36,n29,n31,n34
#SBATCH -o /n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/staging/tenk_xi/tenk_xi.%j.log

set -euo pipefail

WT=/automnt/n17data/cdaley/unions/code/sp_validation.worktrees/repro-paper-ii-astra
PROJ=/n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/Reproductions/UNIONS/pure_eb
SCRIPT=$WT/workflow/scripts/run_2pcf_highres.py
PY=/softs/intelpython/intelpython3-2025.1.0_196/bin/python3.12
OMPI=/softs/openmpi/5.0.5-slurm-CentOS8

export PATH=$OMPI/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:${LD_LIBRARY_PATH:-}
# TreeCorr reads its thread count from SLURM_CPUS_PER_TASK.
export SLURM_CPUS_PER_TASK=12
export OMP_NUM_THREADS=12

NODES=${SLURM_JOB_NUM_NODES:-3}
if [ "${SMOKE:-0}" = "1" ]; then
  NBINS=100; NPATCH=8; PPN=2
  OUT=/n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/staging/tenk_xi/smoke_2node
else
  NBINS=10000; NPATCH=50; PPN=5     # 3 nodes x 5 ranks x 12 threads = 180 cores
  OUT=/n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/staging/tenk_xi/production_out
fi
NTASKS=$(( NODES * PPN ))
mkdir -p "$OUT"

echo "[tenk_xi] $( [ "${SMOKE:-0}" = "1" ] && echo SMOKE || echo PRODUCTION ): nbins=$NBINS npatch=$NPATCH nodes=$NODES ppr=$PPN ntasks=$NTASKS -> $OUT"
echo "[tenk_xi] host mpiexec: $(which mpiexec)"; mpiexec --version | head -1

mpiexec --map-by ppr:${PPN}:node -n ${NTASKS} \
  "$PY" "$SCRIPT" \
  --cat-config $WT/cosmo_val/cat_config.yaml \
  --version    SP_v1.4.6.3_leak_corr \
  --nbins      ${NBINS} --npatch ${NPATCH} \
  --min-sep    0.5 --max-sep 300.0 \
  --out        "$OUT"

echo "[tenk_xi] mpiexec rc=$?  output -> $OUT"
ls -la "$OUT"/*npatch=1.txt 2>/dev/null || echo "[tenk_xi] WARNING: expected txt not found"
