#!/bin/bash
# Usage: run_job.sh <outdir> <nevents> <seed1> <seed2> [particle] [energy] [physlist]
set -e
OUTDIR=$1; NEV=$2; S1=$3; S2=$4; PART=${5:-e-}; EN=${6:-20 GeV}; PL=${7:-FTFP_BERT}
HIDRA=$HOME/HidraSim
source $HIDRA/setup_proof.sh
mkdir -p $OUTDIR && cd $OUTDIR
sed -e "s#^/gps/particle .*#/gps/particle $PART#" -e "s#^/gps/energy .*#/gps/energy $EN#" \
    -e "s#^/random/setSeeds .*#/random/setSeeds $S1 $S2#" -e "s#^/run/beamOn .*#/run/beamOn $NEV#" \
    $HIDRA/DREMTubes_run.mac > run.mac
echo "host=$(hostname) start=$(date)"
$HIDRA/build/DREMTubes -m run.mac -t 1 -pl $PL
echo "end=$(date)"
