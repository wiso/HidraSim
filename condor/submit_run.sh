#!/bin/bash
# Register and submit a sampling-fraction run.
# Usage (on the submit host, proof-11): condor/submit_run.sh <tag> <joblist> [note]
#   env: BUILD=<builddir> (default <repo>/build), PARTICLE=e+, PHYSLIST=FTFP_BERT, PRIO=0 (job priority), MACRO_ARGS='...' (extra make_macro.py args),
#        HIDRA_RUNS=<runs dir>
# Refuses an existing tag. Writes runs/<tag>/run_info.json, updates runs/registry.jsonl and runs/RUNS.md.
set -e
TAG=$1; JOBLIST=$(readlink -f "$2"); NOTE=${3:-}
[ -n "$TAG" ] && [ -f "$JOBLIST" ] || { echo "usage: $0 <tag> <joblist> [note]"; exit 1; }
HERE=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); HIDRA=$(dirname "$HERE")
RUNS=${HIDRA_RUNS:-${OUTBASE:-/storage_tmp/atlas/turra/hidra/sim}}; BUILD=$(readlink -f "${BUILD:-$HIDRA/build}")
PARTICLE=${PARTICLE:-e+}; PHYSLIST=${PHYSLIST:-FTFP_BERT}; PRIO=${PRIO:-0}
[ -e "$RUNS/$TAG" ] && { echo "tag '$TAG' already exists in $RUNS"; exit 1; }
[ -x "$BUILD/DREMTubes" ] || { echo "no executable in $BUILD"; exit 1; }
export HIDRA_DIR=$HIDRA HIDRA_RUNS=$RUNS HIDRA_EXE=$BUILD/DREMTubes
ENVF=$BUILD/hidra_env.sh; [ -f "$ENVF" ] || ENVF=$HIDRA/setup_proof.sh
python3 "$HERE/registry.py" new "$TAG" "$JOBLIST" --note "$NOTE" --particle "$PARTICLE" --physlist "$PHYSLIST" --setup "$ENVF" --macro-args "${MACRO_ARGS:-}"
cd "$HERE"
OUT=$(condor_submit sampling.sub hidra="$HIDRA" runsdir="$RUNS" build="$BUILD" tag="$TAG" \
      joblist="$RUNS/$TAG/joblist.txt" particle="$PARTICLE" physlist="$PHYSLIST" prio="$PRIO" macroargs="${MACRO_ARGS:-}")
echo "$OUT"
CID=$(echo "$OUT" | sed -nE 's/.*submitted to cluster ([0-9]+).*/\1/p')
python3 "$HERE/registry.py" set "$TAG" cluster_id="$CID" status='"running"' condor_priority="$PRIO" submitted="\"$(date -Iseconds)\""
echo "registered $TAG (cluster $CID) in $RUNS/RUNS.md"
