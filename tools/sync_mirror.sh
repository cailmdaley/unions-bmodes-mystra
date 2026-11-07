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

# Copy the report's artifacts from lc's qualified-output layout
# (results/<universe>/<a>/<b>/<file>) into MySTRA's owner layout
# (analyses/<a>/analyses/<b>/results/<universe>/<file>); root outputs stay at
# results/<universe>/<file>. Every file under the size cap is copied,
# which covers the report's figures and values and the inventory's metrics and
# evidence. Inside the sweeps only json, csv and png are copied, so
# per-realisation data stay in the private repository, as do the large
# integration tables. A file whose annex content is absent stays pending.
max_bytes=${MIRROR_MAX_BYTES:-2000000}
copy_artifact() {
  local input=$1 output=$2
  if [[ ! -e "$input" ]] || head -c 15 "$input" | grep -qFx '/annex/objects/'; then
    echo "Pending: ${input#"$source_results"/}" >&2
    return
  fi
  if (( $(stat -L -c %s "$input") > max_bytes )); then
    return
  fi
  mkdir -p "$(dirname "$output")"
  rm -f "$output"
  cp -L -- "$input" "$output"
}
for universe in baseline submitted; do
  while IFS= read -r -d '' input; do
    rel=${input#"$source_results/$universe"/}
    dir=$(dirname "$rel")
    if [[ "$dir" == . ]]; then
      copy_artifact "$input" "$mirror/results/$universe/$rel"
      continue
    fi
    owner="analyses/${dir//\//\/analyses\/}"
    if compgen -G "$project/$owner/*.sweep.yaml" >/dev/null && [[ ! "$rel" =~ \.(json|csv|png)$ ]]; then
      continue
    fi
    copy_artifact "$input" "$mirror/$owner/results/$universe/$(basename "$rel")"
  done < <(find -L "$source_results/$universe" -type f ! -name '.*' ! -path '*/variants/*' -print0)
done

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
