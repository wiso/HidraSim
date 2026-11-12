# Environment of build/ (as used by runs sampling_v1, sampling_v2 and earlier).
# NOTE: this build is linked against Geant4 11.2.2 from the LCG_106b view (see CMakeCache Geant4_DIR);
# the second line only adds the 11.3.1 geant4-config/scripts to PATH and does NOT change the
# Geant4 libraries or data sets used by this binary (they stay 11.2.2 / LCG_106b).
# Kept identical to the old setup_proof.sh so that pending jobs keep the same environment.
source /cvmfs/sft.cern.ch/lcg/views/LCG_106b/x86_64-el9-gcc11-opt/setup.sh
source /cvmfs/geant4.cern.ch/geant4/11.3.p01/x86_64-el9-gcc11-optdeb-MT/bin/geant4.sh
