#!/bin/bash
# Merge a finished run per energy and close its registry entry.
# Usage: condor/close_run.sh <tag> [note]     (env: HIDRA_RUNS, NO_MERGE=1 to skip hadd)
# Merged files: runs/<tag>/merged/sf_E<E>.root (only successful job outputs).
set -e
TAG=$1; NOTE=${2:-}
HERE=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); HIDRA=$(dirname "$HERE")
RUNS=${HIDRA_RUNS:-${OUTBASE:-/storage_tmp/atlas/turra/hidra/sim}}; D=$RUNS/$TAG
[ -f "$D/run_info.json" ] || { echo "unknown tag $TAG"; exit 1; }
export HIDRA_DIR=$HIDRA HIDRA_RUNS=$RUNS
CID=$(python3 -c "import json;print(json.load(open('$D/run_info.json')).get('cluster_id') or '')")
if [ -n "$CID" ] && [ -n "$(condor_q "$CID" -af ProcId 2>/dev/null)" ]; then echo "cluster $CID still has jobs in the queue"; exit 1; fi
MERGED=()
if [ -z "$NO_MERGE" ]; then
  ENVF=$(python3 -c "import json;print(json.load(open('$D/run_info.json'))['environment'].get('setup_file') or '')")
  [ -n "$ENVF" ] && source "$ENVF" >/dev/null 2>&1
  mkdir -p "$D/merged"
  for E in $(awk '{print $1}' "$D/joblist.txt" | sort -un); do
    shopt -s nullglob; fl=("$D"/out/sf_E${E}_*.root "$D"/sf_E${E}_*.root); shopt -u nullglob
    [ ${#fl[@]} -gt 0 ] || continue
    hadd -f "$D/merged/sf_E$E.root" "${fl[@]}" > /dev/null && MERGED+=("merged/sf_E$E.root")
  done
fi
python3 "$HERE/registry.py" close "$TAG" --merged "${MERGED[@]}" ${NOTE:+--note "$NOTE"}
