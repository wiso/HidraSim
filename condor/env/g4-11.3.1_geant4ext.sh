# Environment of build-g4-11.3.1/: Geant4 11.3.1 (geant4.cern.ch 11.3.p01, el9 gcc11 optdeb-MT),
# with the externals view it was built against (LCG_106a_geant4ext20241128: ROOT, CLHEP, XercesC, Qt, gcc11)
# and the Geant4 11.3.1 data sets (geant4.sh does not set them).
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a_geant4ext20241128/x86_64-el9-gcc11-opt/setup.sh
source /cvmfs/geant4.cern.ch/geant4/11.3.p01/x86_64-el9-gcc11-optdeb-MT/bin/geant4.sh
export G4NEUTRONHPDATA=/cvmfs/geant4.cern.ch/share/data/G4NDL4.7.1
export G4LEDATA=/cvmfs/geant4.cern.ch/share/data/G4EMLOW8.6.1
export G4LEVELGAMMADATA=/cvmfs/geant4.cern.ch/share/data/PhotonEvaporation6.1
export G4RADIOACTIVEDATA=/cvmfs/geant4.cern.ch/share/data/RadioactiveDecay6.1.2
export G4PARTICLEXSDATA=/cvmfs/geant4.cern.ch/share/data/G4PARTICLEXS4.1
export G4PIIDATA=/cvmfs/geant4.cern.ch/share/data/G4PII1.3
export G4REALSURFACEDATA=/cvmfs/geant4.cern.ch/share/data/RealSurface2.2
export G4SAIDXSDATA=/cvmfs/geant4.cern.ch/share/data/G4SAIDDATA2.0
export G4ABLADATA=/cvmfs/geant4.cern.ch/share/data/G4ABLA3.3
export G4INCLDATA=/cvmfs/geant4.cern.ch/share/data/G4INCL1.2
export G4ENSDFSTATEDATA=/cvmfs/geant4.cern.ch/share/data/G4ENSDFSTATE3.0
export G4CHANNELINGDATA=/cvmfs/geant4.cern.ch/share/data/G4CHANNELING1.0
export G4PARTICLEHPDATA=/cvmfs/geant4.cern.ch/share/data/G4TENDL1.4
export G4NUDEXLIBDATA=/cvmfs/geant4.cern.ch/share/data/G4NUDEXLIB1.0
export G4URRPTDATA=/cvmfs/geant4.cern.ch/share/data/G4URRPT1.1
