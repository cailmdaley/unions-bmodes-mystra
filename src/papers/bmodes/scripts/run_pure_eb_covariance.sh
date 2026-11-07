#!/usr/bin/env bash
# Pure E/B covariance driver.
#
# --covariance-method monte_carlo runs the 20 independent MC chunks in parallel
# (each a fresh python process, deterministic seed 42+chunk_id) throttled to
# the node core count, then gathers them. --covariance-method exact_linear skips
# the chunks: the gather step builds the exact M C Mᵀ from --cov-integration
# directly. Only the bin centres and pair counts of the two ξ± files enter.
#
# Usage:
#   run_pure_eb_covariance.sh --version SP_v1.4.6.3_leak_corr --blind A \
#     --cat-config <cat_config.yaml> \
#     --xi-reporting <xi 20-bin .txt> --xi-integration <xi 1000-bin .txt> \
#     --cov-integration <cov ..._processed.txt> \
#     --integration-method {adaptive_dqags|fixed_operator} \
#     --integration-limits {floor_1arcmin|fine_grid_extent|reporting_extent} \
#     --rebinning-weights {uniform|npairs|none} \
#     --covariance-method {exact_linear|monte_carlo} \
#     [--geometry {pointwise|averaged}] [--local-from-eval-grid] \
#     --out-file <covariance.npz> \
#     [--n-chunks 20] [--n-samples 2000] [--nproc 16]
set -euo pipefail

# SCRIPTS is this project's vendored scripts; the sp_validation package is
# installed in the recipe's environment.
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VERSION=""; BLIND="A"; CATCONFIG=""; XIREP=""; XIINT=""; COVINT=""; OUTFILE=""
# No default: the pure_eb_integration_method and pure_eb_rebinning_weights
# decisions are the source of truth.
INTMETHOD=""; INTLIMITS=""; REBINWEIGHTS=""; COVMETHOD=""
# Optional: override the geometry the method defaults to. Empty = the default
# mapping.
GEOMETRY=""
# Optional: pointwise fixed_operator only -- read the local 0.5*(xi_p +/- xi_m)
# term from the reporting grid instead of the snapped fine node. Empty = off.
LOCALEVAL=""
# NPROC is how many MC chunk processes this ONE task forks. It must be the
# task's slot, not the allocation. Same knob as the TreeCorr cap in
# run_2pcf.py, set by the runner; --nproc overrides.
NCHUNKS=20; NSAMPLES=""; NPROC="${OMP_NUM_THREADS:-1}"
MINSEP=1.0; MAXSEP=250.0; NBINS=20
MINSEPINT=0.5; MAXSEPINT=300.0; NBINSINT=1000; NPATCH=1

while [ $# -gt 0 ]; do
  case "$1" in
    --version) VERSION="$2"; shift 2;;
    --blind) BLIND="$2"; shift 2;;
    --cat-config) CATCONFIG="$2"; shift 2;;
    --xi-reporting) XIREP="$2"; shift 2;;
    --xi-integration) XIINT="$2"; shift 2;;
    --cov-integration) COVINT="$2"; shift 2;;
    --integration-method) INTMETHOD="$2"; shift 2;;
    --integration-limits) INTLIMITS="$2"; shift 2;;
    --rebinning-weights) REBINWEIGHTS="$2"; shift 2;;
    --covariance-method) COVMETHOD="$2"; shift 2;;
    --geometry) GEOMETRY="$2"; shift 2;;
    --local-from-eval-grid) LOCALEVAL="--local-from-eval-grid"; shift 1;;
    --out-file) OUTFILE="$2"; shift 2;;
    --n-chunks) NCHUNKS="$2"; shift 2;;
    --n-samples) NSAMPLES="$2"; shift 2;;
    --nproc) NPROC="$2"; shift 2;;
    --min-sep) MINSEP="$2"; shift 2;;
    --max-sep) MAXSEP="$2"; shift 2;;
    --nbins) NBINS="$2"; shift 2;;
    --min-sep-int) MINSEPINT="$2"; shift 2;;
    --max-sep-int) MAXSEPINT="$2"; shift 2;;
    --nbins-int) NBINSINT="$2"; shift 2;;
    --npatch) NPATCH="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$INTMETHOD" ] || { echo "--integration-method is required" >&2; exit 2; }
[ -n "$INTLIMITS" ] || { echo "--integration-limits is required" >&2; exit 2; }
[ -n "$REBINWEIGHTS" ] || { echo "--rebinning-weights is required" >&2; exit 2; }
[ -n "$COVMETHOD" ] || { echo "--covariance-method is required" >&2; exit 2; }
[ -n "$NSAMPLES" ] || { echo "--n-samples is required" >&2; exit 2; }
GEOMETRY_FLAG=""
[ -z "$GEOMETRY" ] || GEOMETRY_FLAG="--geometry $GEOMETRY"
LOCALEVAL_FLAG="$LOCALEVAL"

[ -n "$OUTFILE" ] || { echo "--out-file is required" >&2; exit 2; }
CHUNKDIR="$(mktemp -d)"
trap 'rm -rf "$CHUNKDIR"' EXIT

if [ "$COVMETHOD" = monte_carlo ]; then
  echo "[pure_eb] $NCHUNKS chunks, $NSAMPLES samples, nproc=$NPROC, version=$VERSION blind=$BLIND"
  for i in $(seq 0 $((NCHUNKS-1))); do
    (
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        python "$SCRIPTS/precompute_pure_eb_chunk.py" \
        --chunk-id "$i" --n-chunks "$NCHUNKS" --n-samples "$NSAMPLES" \
        --version "$VERSION" --blind "$BLIND" --cat-config "$CATCONFIG" \
        --xi-reporting "$XIREP" --xi-integration "$XIINT" --cov-integration "$COVINT" \
        --min-sep "$MINSEP" --max-sep "$MAXSEP" --nbins "$NBINS" \
        --min-sep-int "$MINSEPINT" --max-sep-int "$MAXSEPINT" --nbins-int "$NBINSINT" \
        --npatch "$NPATCH" $GEOMETRY_FLAG $LOCALEVAL_FLAG \
        --integration-method "$INTMETHOD" \
        --integration-limits "$INTLIMITS" --rebinning-weights "$REBINWEIGHTS" \
        --out "$CHUNKDIR"
    ) > "$CHUNKDIR/chunk_${i}.log" 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 1; done
  done
  wait

  # Verify all chunks landed
  missing=0
  for i in $(seq 0 $((NCHUNKS-1))); do
    [ -f "$CHUNKDIR/pure_eb_chunk_${i}.npz" ] || { echo "MISSING chunk $i (see $CHUNKDIR/chunk_${i}.log)" >&2; missing=1; }
  done
  [ "$missing" -eq 0 ] || { echo "[pure_eb] chunk failures — aborting gather" >&2; exit 1; }
  echo "[pure_eb] all $NCHUNKS chunks done; gathering"
  COVSOURCE=(--chunks-dir "$CHUNKDIR")
else
  echo "[pure_eb] exact M C M^T covariance, version=$VERSION blind=$BLIND"
  COVSOURCE=(--cov-integration "$COVINT")
fi

python "$SCRIPTS/gather_pure_eb_chunks.py" \
  --xi-reporting "$XIREP" --xi-integration "$XIINT" \
  --covariance-method "$COVMETHOD" "${COVSOURCE[@]}" \
  --min-sep "$MINSEP" --max-sep "$MAXSEP" --nbins "$NBINS" \
  --min-sep-int "$MINSEPINT" --max-sep-int "$MAXSEPINT" --nbins-int "$NBINSINT" \
  $GEOMETRY_FLAG \
  --integration-method "$INTMETHOD" \
  --integration-limits "$INTLIMITS" --rebinning-weights "$REBINWEIGHTS" \
  --out-file "$OUTFILE"
echo "[pure_eb] done -> $OUTFILE"
