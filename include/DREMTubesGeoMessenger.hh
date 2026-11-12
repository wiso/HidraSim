//**************************************************
// \file DREMTubesGeoMessenger.hh
// \brief: Definition of DREMTubesGeomMessenger class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 10 August 2023
//**************************************************

#ifndef DREMTubesGeoMessenger_h
#define DREMTubesGeoMessenger_h 1

//Includers from Geant4
//
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//Includers from project files
//
class DREMTubesDetectorConstruction;

class DREMTubesGeoMessenger final : public G4UImessenger {

    public:
        //Constructor and destructor
        DREMTubesGeoMessenger(DREMTubesDetectorConstruction *DetConstruction);
        ~DREMTubesGeoMessenger();
        
        //Virtual methods from base class
        void SetNewValue(G4UIcommand *command, G4String newValue) override;

    private:
        //Members
        DREMTubesDetectorConstruction *fDetConstruction;
        G4UIdirectory *fMsgrDirectory;
        G4UIcmdWithADoubleAndUnit *fXshiftcmd;
        G4UIcmdWithADoubleAndUnit *fYshiftcmd;
        G4UIcmdWithADoubleAndUnit *fOrzrotcmd;
        G4UIcmdWithADoubleAndUnit *fVerrotcmd;
        G4UIcmdWithAString *fFiberMapcmd;
        G4UIcmdWithAString *fTowerMapcmd;

};

#endif //DREMTubesGeoMessenger_h

//**************************************************
