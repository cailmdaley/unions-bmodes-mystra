#!/usr/bin/env bash
# Stage the committed report and its materialized artifacts in a mirror clone.
set -euo pipefail
if [[ $# != 1 ]]; then
  echo "usage: $0 <mirror-clone>" >&2
  exit 2
fi
project=$(cd "$(dirname "$0")/.." && pwd -P)
repo=$(git -C "$project" rev-parse --show-toplevel)
mirror=$(cd "$1" && pwd -P)
if [[ ! -d "$mirror/.git" || "$mirror" == "$repo" || "$mirror" == "$project" || "$mirror" == "$repo/"* ]]; then
  echo "Refusing to wipe: target must be a separate git clone with a .git directory" >&2
  exit 2
fi
source_results=${LC_RESULTS_ROOT:-/automnt/n17data/cdaley/unions/analyses/shear_2d/bmodes_2d/Reproductions/UNIONS/pure_eb/results}
source_results=$(realpath "$source_results")
echo "Staging $(git -C "$repo" rev-parse --short HEAD) into $mirror"
find "$mirror" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf -- {} +
git -C "$repo" archive HEAD:UNIONS/pure_eb | tar -x -C "$mirror"
find "$mirror" -name '*cov_tmp*' \( -type f -o -type l \) -delete

# src is lc's qualified-output layout; dst is the MySTRA SDK's owner layout.
# Only exact declared formats are copied: missing PNGs must remain pending.
while IFS='|' read -r src dst; do
  [[ -z "$src" ]] && continue
  input="$source_results/$src"
  output="$mirror/$dst"
  if [[ ! -f "$input" ]]; then
    echo "Pending: $src" >&2
    continue
  fi
  if head -c 15 "$input" | grep -qFx '/annex/objects/'; then
    echo "Pending (annex pointer): $src" >&2
    continue
  fi
  mkdir -p "$(dirname "$output")"
  rm -f "$output"
  cp -L -- "$input" "$output"
done <<'ARTIFACTS'
baseline/null_tests/paper_claims.csv|analyses/null_tests/results/baseline/paper_claims.csv
baseline/null_tests/paper_claims_comparison.csv|analyses/null_tests/results/baseline/paper_claims_comparison.csv
submitted/null_tests/paper_claims.csv|results/submitted/null_tests/paper_claims.csv
baseline/pte_table_appendix.csv|results/baseline/pte_table_appendix.csv
baseline/pure_eb/pure_eb_data_vector_figure.png|analyses/pure_eb/results/baseline/pure_eb_data_vector_figure.png
baseline/pure_eb/pure_eb_covariance_figure.png|analyses/pure_eb/results/baseline/pure_eb_covariance_figure.png
baseline/pure_eb/pure_eb_version_comparison_figure.png|analyses/pure_eb/results/baseline/pure_eb_version_comparison_figure.png
baseline/cosebis/fig_cosebis_data_vector.png|analyses/cosebis/results/baseline/fig_cosebis_data_vector.png
baseline/cosebis/cosebis_version_comparison.png|analyses/cosebis/results/baseline/cosebis_version_comparison.png
baseline/cl_bb/cl_data_vector_figure.png|analyses/cl_bb/results/baseline/cl_data_vector_figure.png
baseline/cl_bb/cl_version_comparison_figure.png|analyses/cl_bb/results/baseline/cl_version_comparison_figure.png
baseline/glass_mocks/mock_validation_figure.png|analyses/glass_mocks/results/baseline/mock_validation_figure.png
baseline/glass_mocks/mock_min_pte_figure.png|analyses/glass_mocks/results/baseline/mock_min_pte_figure.png
baseline/fig_pte_heatmaps.png|results/baseline/fig_pte_heatmaps.png
baseline/fig_cl_pte_heatmap.png|results/baseline/fig_cl_pte_heatmap.png
baseline/fig_harmonic_config_cosebis_full.png|results/baseline/fig_harmonic_config_cosebis_full.png
ARTIFACTS

# Remove unfetched annex pointers and annex symlinks, without following links.
python3 - "$mirror" <<'PY'
import os, sys
from pathlib import Path
root = Path(sys.argv[1])
for directory, dirs, files in os.walk(root, followlinks=False):
    dirs[:] = [d for d in dirs if d != '.git']
    for name in files:
        path = Path(directory, name)
        if path.is_symlink():
            if '/annex/objects/' in os.readlink(path):
                path.unlink()
        elif path.is_file():
            with path.open('rb') as f:
                if f.read(15) == b'/annex/objects/':
                    path.unlink()
PY
mkdir -p "$mirror/results"
printf '%s\n' 'The full artifacts and their manifests live in the private reproduction repository.' > "$mirror/results/README.md"
echo "Mirror size (excluding .git):"
du -sh --exclude=.git "$mirror"
