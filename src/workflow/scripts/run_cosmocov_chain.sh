#!/usr/bin/env bash
# CosmoCov covariance chain.
#
# Generates the CosmoCov .ini, runs the 3 shear-shear blocks (++,--,+-), which
# are independent, in parallel, concatenates them and post-processes the
# result with cosmocov_process. The .ini generation and cosmocov_process steps
# run with the ambient python (sp_validation is installed in the project
# environment).
#
# The CosmoCov binary is dynamically linked against libgsl.so.23,
# libgslcblas and libfftw3. The container carries the latter two
# but ships gsl.so.27, so the one missing soname is supplied from the
# host via --libgsl and linked into a scratch directory that is on
# LD_LIBRARY_PATH *only* for the binary's own invocations — letting
# gsl 23 shadow the container stack for the python steps would put
# pyccl and friends on a library they were not built
# against.
#
# The whole chain runs in a private scratch directory: no CosmoCov temporaries
# (cov_tmp*/order_cov_tmp*), .ini or .log ever land next to the output. Exactly
# one product, selected by --emit, is copied to --out-file:
#
#   blocks       concatenated raw ++/--/+- blocks   (covariance.txt)
#   processed    processed G+NG matrix              (covariance_processed.txt)
#   processed_g  Gaussian-only extract              (covariance_processed_g.txt)
#   plot         correlation-structure QA figure    (300-dpi PNG)
#
# Usage:
#   run_cosmocov_chain.sh --version SP_v1.4.6.3_leak_corr --blind blind_a \
#     --min-sep 0.5 --max-sep 300.0 --nbins 1000 --gaussian g \
#     --planck18-json <cosmology_snapshot.json> \
#     --mask-cls <mask_cls...norm.txt> \
#     [--nz <n(z) file; default the per-version SOM n(z) by path convention>] \
#     --cosmocov <cov binary> --libgsl <libgsl.so.23.x> \
#     --emit processed --out-file <path>
set -euo pipefail

# VENDOR_ROOT is this project's vendored script tree (src/).
VENDOR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COSMOCOV=/n23data1/n06data/lgoh/scratch/UNIONS/CosmoCov/covs/cov
LIBGSL=/lib64/libgsl.so.23.1.0

VERSION=""; BLIND=""; MINSEP=""; MAXSEP=""; NBINS=""; GAUSSIAN=""
# A supplied spectrum applies to SSC unless explicitly disabled.
PLANCK18=""; MASKCLS=""; MASKAPPLICATION="masked"; NZ=""; EMIT=""; OUTFILE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --version) VERSION="$2"; shift 2;;
    --blind) BLIND="$2"; shift 2;;
    --min-sep) MINSEP="$2"; shift 2;;
    --max-sep) MAXSEP="$2"; shift 2;;
    --nbins) NBINS="$2"; shift 2;;
    --gaussian) GAUSSIAN="$2"; shift 2;;
    --planck18-json) PLANCK18="$2"; shift 2;;
    --mask-cls) MASKCLS="$2"; shift 2;;
    --mask-application) MASKAPPLICATION="$2"; shift 2;;
    --nz) NZ="$2"; shift 2;;
    --cosmocov) COSMOCOV="$2"; shift 2;;
    --libgsl) LIBGSL="$2"; shift 2;;
    --emit) EMIT="$2"; shift 2;;
    --out-file) OUTFILE="$2"; shift 2;;
    -h|--help) sed -n '2,34p' "${BASH_SOURCE[0]}"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

# ASTRA option IDs and the single-letter spelling name the same terms.
case "$GAUSSIAN" in
  gaussian_only) GAUSSIAN=g;;
esac
case "$MASKAPPLICATION" in
  masked) ;;
  unmasked) MASKCLS="";;
  *) echo "--mask-application must be masked or unmasked" >&2; exit 2;;
esac
# The nside-4096 spectra are too short for Gaussian w_mask() on either grid.
# They affect only survey_variance() in SSC; do not pass an unused file for G.
if [ "$GAUSSIAN" = g ]; then MASKCLS=""; fi

case "$EMIT" in
  blocks|processed|processed_g|plot) ;;
  *) echo "--emit must be one of blocks|processed|processed_g|plot" >&2; exit 2;;
esac
[ -n "$OUTFILE" ] || { echo "--out-file is required" >&2; exit 2; }

# Absolutize the output path before any `cd`: the CosmoCov binary writes its
# block files into cwd, so the chain cds into the scratch dir below.
OUTDIR="$(cd "$(dirname "$OUTFILE")" && pwd)"
OUTFILE="$OUTDIR/$(basename "$OUTFILE")"

# Private scratch: everything CosmoCov emits (cov_tmp*, order_cov_tmp*, .ini,
# .log) stays here and is removed on exit.
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/cosmocov.XXXXXXXX")"
trap 'rm -rf "$SCRATCH"' EXIT
INI="$SCRATCH/covariance.ini"

# The binary reads the mask and n(z) files the .ini names relative to its cwd
# (the scratch dir) and truncates long paths, so both are copied in and named
# by basename.
if [ -n "$MASKCLS" ]; then cp "$MASKCLS" "$SCRATCH/mask_cls.txt"; MASKCLS=mask_cls.txt; fi
if [ -n "$NZ" ]; then cp "$NZ" "$SCRATCH/nz.txt"; NZ=nz.txt; fi

echo "[cosmocov] generating .ini"
python "$VENDOR_ROOT/workflow/scripts/generate_cosmocov_ini.py" \
  --version "$VERSION" --blind "$BLIND" \
  --planck18-json "$PLANCK18" \
  --min-sep "$MINSEP" --max-sep "$MAXSEP" --nbins "$NBINS" --gaussian "$GAUSSIAN" \
  --mask-cls "$MASKCLS" ${NZ:+--nz "$NZ"} --out-ini "$INI"

# The binary resolves libgsl.so.23 out of here; nothing else sees it.
HOSTLIBS="$SCRATCH/hostlibs"
mkdir -p "$HOSTLIBS"
ln -sf "$LIBGSL" "$HOSTLIBS/libgsl.so.23"

echo "[cosmocov] running 3 blocks (parallel)"
cd "$SCRATCH"
# One CosmoCov invocation per block: ++ is 1, -- is 2, +- is 3
for idx in 1 2 3; do
  (
    LD_LIBRARY_PATH="$HOSTLIBS${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
      "$COSMOCOV" "$idx" "$INI" > "$SCRATCH/cosmocov_block_${idx}.log" 2>&1
  ) &
done
wait

# Concatenate blocks in order (++, --, +-)
CAT="$SCRATCH/covariance.txt"
: > "$CAT"
for pm_idx in "++:1" "--:2" "+-:3"; do
  pm="${pm_idx%%:*}"; idx="${pm_idx##*:}"
  blk="$SCRATCH/cov_tmp_ssss_${pm}_cov_Ntheta${NBINS}_Ntomo1_${idx}"
  [ -f "$blk" ] || { echo "MISSING block $blk (see $SCRATCH/cosmocov_block_${idx}.log)" >&2; exit 1; }
  cat "$blk" >> "$CAT"
done
echo "[cosmocov] concatenated -> $CAT"

if [ "$EMIT" = "blocks" ]; then
  cp "$CAT" "$OUTFILE"
  echo "[cosmocov] done -> $OUTFILE"
  exit 0
fi

echo "[cosmocov] processing (positive-definite check, G/G+NG extract, QA plot)"
# Third argument is the explicit plot path: a .png is rendered at 300 dpi.
python "$VENDOR_ROOT/cosmo_inference/scripts/cosmocov_process.py" \
  "$CAT" "$SCRATCH/covariance_processed" "$SCRATCH/covariance_processed_plot.png"

case "$EMIT" in
  processed)   PRODUCT="$SCRATCH/covariance_processed.txt";;
  processed_g) PRODUCT="$SCRATCH/covariance_processed_g.txt";;
  plot)        PRODUCT="$SCRATCH/covariance_processed_plot.png";;
esac
cp "$PRODUCT" "$OUTFILE"
echo "[cosmocov] done -> $OUTFILE"
