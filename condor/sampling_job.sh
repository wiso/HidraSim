#!/bin/bash
# Sampling-fraction study job (one HTCondor job, 1 thread, runs in place on the shared FS).
# Usage: sampling_job.sh <rundir> <energy_GeV> <iE> <job> <nevents> [builddir] [particle] [physlist] [make_macro args...]
#   seed1 = 1000*iE + job + 1 (seed2 = 7*seed1+13, see sampling/make_macro.py)
#   environment: <builddir>/hidra_env.sh (falls back to $HIDRA_SETUP or <repo>/setup_proof.sh)
#   output: <rundir>/sf_E<energy>_<job>.root, job record <rundir>/jobinfo/E<energy>_<job>.json
set -e
RUNDIR=$1; EN=$2; IE=$3; JOB=$4; NEV=$5
HIDRA=${HIDRA_DIR:-$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)}
BUILD=${6:-$HIDRA/build}; PART=${7:-e+}; PL=${8:-FTFP_BERT}
shift $(( $# < 8 ? $# : 8 )); MACARGS="$*"
EXE=$BUILD/DREMTubes
ENVF=$BUILD/hidra_env.sh
[ -f "$ENVF" ] || ENVF=${HIDRA_SETUP:-$HIDRA/setup_proof.sh}
source "$ENVF"
SEED=$((1000*IE + JOB + 1))
WORK=$RUNDIR/work/E${EN}_${JOB}
mkdir -p "$WORK" "$RUNDIR/jobinfo" && cd "$WORK"
python3 "$HIDRA/sampling/make_macro.py" --energy "$EN" --nevents "$NEV" --seed "$SEED" --particle "$PART" $MACARGS > job.mac
T0=$(date +%s)
echo "host=$(hostname) E=$EN job=$JOB seed=$SEED build=$BUILD env=$ENVF start=$(date)"
set +e
"$EXE" -m job.mac -t 1 -pl "$PL" 2>&1 | tee run.log
RC=${PIPESTATUS[0]}
set -e
T1=$(date +%s)
G4=$(grep -m1 -oE 'geant4-[0-9]+-[0-9]+(-patch-[0-9]+)?' run.log || true)
G4LIB=$(ldd "$EXE" | awk '/libG4run/{print $3}')
OS=$(grep PRETTY_NAME /etc/os-release | cut -d= -f2 | tr -d '"')
DATA=$(env | grep -E '^G4[A-Z]*DATA=' | sed -E 's#^([A-Z0-9]+)=.*/([^/]+)$#"\1":"\2"#' | paste -sd, -)
[ $RC -eq 0 ] && [ -f DREMTubesout_Run0.root ] && mv DREMTubesout_Run0.root "$RUNDIR/sf_E${EN}_${JOB}.root"
cat > "$RUNDIR/jobinfo/E${EN}_${JOB}.json" <<J
{"energy":"$EN","job":$JOB,"ie":$IE,"seed1":$SEED,"nev":$NEV,"particle":"$PART","physics_list":"$PL","macro_args":"$MACARGS",
 "host":"$(hostname)","os":"$OS","start":$T0,"end":$T1,"elapsed_s":$((T1-T0)),"exit_code":$RC,
 "geant4":"$G4","geant4_lib":"$G4LIB","executable":"$EXE","executable_md5":"$(md5sum "$EXE" | cut -d' ' -f1)",
 "env_file":"$ENVF","g4_datasets":{$DATA},"condor_job":"${_CONDOR_JOB_AD:+$(awk -F' = ' '/^GlobalJobId/{gsub(/"/,"",$2);print $2}' "$_CONDOR_JOB_AD")}"}
J
echo "end=$(date) elapsed_s=$((T1-T0)) nev=$NEV exit=$RC geant4=$G4"
exit $RC
