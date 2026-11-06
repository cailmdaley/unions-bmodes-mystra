#!/usr/bin/env bash
# Pure E/B PTE-matrix sweep over non-fiducial versions (Design B, Tier-2).
#
# Loops the non-fiducial version list (from catalog_versions.py --nonfiducial) and runs
# calculate_pure_eb_ptes.py once per version — the same χ² PTE-matrix compute the
# fiducial pure_eb_pte_per_cut recipe calls. Each version's gathered semi-analytic
# NPZ is read from the pure_eb_sweep output dir and its Gaussian integration
# covariance from the cov_sweep output dir (both by absolute path — lc does not
# wire cross-output deps, so the driver reads upstream sweeps directly; run
# pure_eb_sweep + cov_sweep before this). Each version emits the canonical
# ``{ver}_{blind}_pure_eb_ptes.npz`` — the exact name config_space_pte_matrices.py
# reconstructs from --pte-intermediate-dir — straight into --out. Serial over
# versions; each version is fast (206-pair grid, ~seconds).
#
# Usage:
#   run_pure_eb_ptes_sweep.sh --cat-config <cat_config.yaml> \
#     --pure-eb-sweep-dir <pure_eb/.../pure_eb_sweep> \
#     --cov-sweep-dir <covariance/.../cov_sweep> \
#     --out <output_dir> --blind A --versions "v1 v2 ..."
set -euo pipefail

# WSCRIPTS/PSCRIPTS are this project's vendored scripts; the sp_validation
# package is installed in the environment lc runs the recipe in.
VENDOR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
WSCRIPTS=$VENDOR_ROOT/workflow/scripts
PSCRIPTS=$VENDOR_ROOT/papers/bmodes/scripts

CATCONFIG=""; PUREEBSWEEP=""; COVSWEEP=""; OUT=""; BLIND=""; VERSIONS=""
while [ $# -gt 0 ]; do
  case "$1" in
    --cat-config) CATCONFIG="$2"; shift 2;;
    --pure-eb-sweep-dir) PUREEBSWEEP="$2"; shift 2;;
    --cov-sweep-dir) COVSWEEP="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --blind) BLIND="$2"; shift 2;;
    --versions) VERSIONS="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$VERSIONS" ] || { echo "--versions is required" >&2; exit 2; }

mkdir -p "$OUT"

for ver in $VERSIONS; do
  pureeb="$PUREEBSWEEP/${ver}_${BLIND}_pure_eb_semianalytic.npz"
  covbase="covariance_${ver}_${BLIND}_integration"
  covint="$COVSWEEP/$covbase/${covbase}_processed.txt"
  for f in "$pureeb" "$covint"; do
    [ -f "$f" ] || { echo "MISSING upstream input for $ver: $f" >&2; exit 1; }
  done
  echo "[pure_eb_ptes_sweep] $ver"
  python "$PSCRIPTS/calculate_pure_eb_ptes.py" \
    --version "$ver" --blind "$BLIND" \
    --pure-eb-data "$pureeb" --cov-integration "$covint" \
    --npatch 1 --out "$OUT"
  echo "[pure_eb_ptes_sweep] $ver -> $OUT/${ver}_${BLIND}_pure_eb_ptes.npz"
done
echo "[pure_eb_ptes_sweep] done -> $OUT"
