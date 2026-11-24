#!/usr/bin/env bash
# CosmoCov covariance chain.
#
# Generates the CosmoCov .ini, runs the 3 shear-shear blocks (++,--,+-), which
# are independent, in parallel, concatenates them and post-processes the
# result with cosmocov_process. The .ini generation and cosmocov_process steps
# run with the ambient python (sp_validation is installed in the project
# environment).
#
# The CosmoCov binary is built into the project image from the public
# CosmoLike/CosmoCov code (pyproject.toml, [tool.lightcone.image]).
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
#     [--nz <n(z) file; default the blind's decision_values.NZ_PATH>] \
#     [--reference-covariance <option setting the survey block; default the
#      version's own survey properties>] \
#     --emit processed --out-file <path>
set -euo pipefail

# VENDOR_ROOT is this project's vendored script tree (src/).
VENDOR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COSMOCOV=/opt/CosmoCov/covs/cov

VERSION=""; BLIND=""; MINSEP=""; MAXSEP=""; NBINS=""; GAUSSIAN=""
PLANCK18=""; NZ=""; REFERENCE=""; EMIT=""; OUTFILE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --version) VERSION="$2"; shift 2;;
    --blind) BLIND="$2"; shift 2;;
    --min-sep) MINSEP="$2"; shift 2;;
    --max-sep) MAXSEP="$2"; shift 2;;
    --nbins) NBINS="$2"; shift 2;;
    --gaussian) GAUSSIAN="$2"; shift 2;;
    --planck18-json) PLANCK18="$2"; shift 2;;
    --nz) NZ="$2"; shift 2;;
    --reference-covariance) REFERENCE="$2"; shift 2;;
    --cosmocov) COSMOCOV="$2"; shift 2;;
    --emit) EMIT="$2"; shift 2;;
    --out-file) OUTFILE="$2"; shift 2;;
    -h|--help) sed -n '2,29p' "${BASH_SOURCE[0]}"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

# ASTRA option IDs and the single-letter spelling name the same terms.
case "$GAUSSIAN" in
  gaussian_only) GAUSSIAN=g;;
esac

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

# The binary reads the n(z) file the .ini names relative to its cwd (the scratch
# dir) and truncates long paths, so the n(z) is copied in as nz.txt.
[ -n "$NZ" ] || NZ="$(python "$VENDOR_ROOT/workflow/scripts/decision_values.py" nz_file "$BLIND")"
cp "$NZ" "$SCRATCH/nz.txt"; NZ=nz.txt

echo "[cosmocov] generating .ini"
python "$VENDOR_ROOT/workflow/scripts/generate_cosmocov_ini.py" \
  --version "$VERSION" --blind "$BLIND" \
  --planck18-json "$PLANCK18" \
  --min-sep "$MINSEP" --max-sep "$MAXSEP" --nbins "$NBINS" --gaussian "$GAUSSIAN" \
  --nz "$NZ" \
  ${REFERENCE:+--reference-covariance "$REFERENCE"} --out-ini "$INI"

echo "[cosmocov] running 3 blocks (parallel)"
cd "$SCRATCH"
# One CosmoCov invocation per block: ++ is 1, -- is 2, +- is 3
for idx in 1 2 3; do
  "$COSMOCOV" "$idx" "$INI" > "$SCRATCH/cosmocov_block_${idx}.log" 2>&1 &
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
