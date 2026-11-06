#!/usr/bin/env bash
# Pure E/B semi-analytic covariance sweep over non-fiducial versions (Design B).
#
# Loops the non-fiducial version list (from catalog_versions.py --nonfiducial) and runs
# the same run_pure_eb_semianalytic.sh scatter-gather the fiducial
# pure_eb_semianalytic_data recipe calls, once per version. Each version's
# reporting+integration ξ± is read from the xi_sweep output dir and its Gaussian
# integration covariance from the cov_sweep output dir (both by absolute path —
# lc does not wire cross-output deps, so the driver reads upstream sweeps
# directly; run xi_sweep + cov_sweep before this). The gathered
# ``{ver}_{blind}_pure_eb_semianalytic.npz`` — already the canonical name
# pure_eb_version_comparison (--results-dir) reads — is moved into --out; the
# per-version MC chunks stage in a scratch subdir that is removed after gather.
#
# Usage:
#   run_pure_eb_sweep.sh --cat-config <cat_config.yaml> \
#     --xi-sweep-dir <two_point/.../xi_sweep> \
#     --cov-sweep-dir <covariance/.../cov_sweep> \
#     --integration-method {adaptive_dqags|fixed_operator} \
#     --integration-limits {floor_1arcmin|fine_grid_extent|reporting_extent} \
#     --rebinning-weights {uniform|npairs|none} \
#     --covariance-method {exact_linear|monte_carlo} \
#     --out <output_dir> --blind A --versions "v1 v2 ..."
set -euo pipefail

# WSCRIPTS/PSCRIPTS are this project's vendored scripts; the sp_validation
# package is installed in the environment lc runs the recipe in.
VENDOR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
WSCRIPTS=$VENDOR_ROOT/workflow/scripts
PSCRIPTS=$VENDOR_ROOT/papers/bmodes/scripts

CATCONFIG=""; XISWEEP=""; COVSWEEP=""; OUT=""; BLIND=""; VERSIONS=""
# No default: the pure_eb_integration_method and pure_eb_rebinning_weights
# decisions are the source of truth.
INTMETHOD=""; INTLIMITS=""; REBINWEIGHTS=""; COVMETHOD=""; NSAMPLES=""
while [ $# -gt 0 ]; do
  case "$1" in
    --cat-config) CATCONFIG="$2"; shift 2;;
    --xi-sweep-dir) XISWEEP="$2"; shift 2;;
    --cov-sweep-dir) COVSWEEP="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --blind) BLIND="$2"; shift 2;;
    --versions) VERSIONS="$2"; shift 2;;
    --integration-method) INTMETHOD="$2"; shift 2;;
    --integration-limits) INTLIMITS="$2"; shift 2;;
    --rebinning-weights) REBINWEIGHTS="$2"; shift 2;;
    --covariance-method) COVMETHOD="$2"; shift 2;;
    --n-samples) NSAMPLES="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$INTMETHOD" ] || { echo "--integration-method is required" >&2; exit 2; }
[ -n "$INTLIMITS" ] || { echo "--integration-limits is required" >&2; exit 2; }
[ -n "$REBINWEIGHTS" ] || { echo "--rebinning-weights is required" >&2; exit 2; }
[ -n "$COVMETHOD" ] || { echo "--covariance-method is required" >&2; exit 2; }
[ -n "$NSAMPLES" ] || { echo "--n-samples is required" >&2; exit 2; }

[ -n "$VERSIONS" ] || { echo "--versions is required" >&2; exit 2; }

mkdir -p "$OUT"

for ver in $VERSIONS; do
  xirep="$XISWEEP/${ver}_xi_minsep=1.0_maxsep=250.0_nbins=20_npatch=1.txt"
  xiint="$XISWEEP/${ver}_xi_minsep=0.5_maxsep=300.0_nbins=1000_npatch=1.txt"
  covbase="covariance_${ver}_${BLIND}_integration"
  covint="$COVSWEEP/$covbase/${covbase}_processed.txt"
  for f in "$xirep" "$xiint" "$covint"; do
    [ -f "$f" ] || { echo "MISSING upstream input for $ver: $f" >&2; exit 1; }
  done
  echo "[pure_eb_sweep] $ver"
  stage="$OUT/_stage_${ver}"
  bash "$PSCRIPTS/run_pure_eb_semianalytic.sh" \
    --version "$ver" --blind "$BLIND" --cat-config "$CATCONFIG" \
    --xi-reporting "$xirep" --xi-integration "$xiint" --cov-integration "$covint" \
    --integration-method "$INTMETHOD" \
    --integration-limits "$INTLIMITS" \
    --rebinning-weights "$REBINWEIGHTS" \
    --covariance-method "$COVMETHOD" --n-samples "$NSAMPLES" \
    --out "$stage"
  mv "$stage/${ver}_${BLIND}_pure_eb_semianalytic.npz" \
     "$OUT/${ver}_${BLIND}_pure_eb_semianalytic.npz"
  rm -rf "$stage"
  echo "[pure_eb_sweep] $ver -> $OUT/${ver}_${BLIND}_pure_eb_semianalytic.npz"
done
echo "[pure_eb_sweep] done -> $OUT"
