//**************************************************
// \file DREMTubesDetectorConstruction.hh
// \brief: Definition of DREMTubesDetectorConstruction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multplie times
//
#ifndef DREMTubesDetectorConstruction_h
#define DREMTubesDetectorConstruction_h 1

//Includers from Geant4
//
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

//Includers from C++
//
#include <vector>
#include <string>

//Include geometrical parameters
#include "DREMTubesGeoPar.hh"
#include "DREMTubesGeoMessenger.hh"
//Forward declaration
//
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class DREMTubesDetectorConstruction : public G4VUserDetectorConstruction {
  
    public:
        //Constructor
        //
        DREMTubesDetectorConstruction();
        //De-constructor
        //
        virtual ~DREMTubesDetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct();

        G4LogicalVolume* constructscinfiber(double tolerance,
                                            G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* ScinMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CherMaterial);
    
        G4LogicalVolume* constructcherfiber(double tolerance, 
                                            G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* CherMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CladCherMaterial);
        //Getters
    	//
        //const G4VPhysicalVolume* GetLeakCntPV() const;
        const G4VPhysicalVolume* GetLeakCntlPV() const;
        const G4VPhysicalVolume* GetLeakCntdPV() const;
        const G4VPhysicalVolume* GetWorldPV() const;

        //Other methods
        //
        G4int GetTowerID( const G4int& cpno ) const;
        G4int GetSiPMID(const G4int& cpno ) const; 
        G4int GetSiPMTower(const G4int& town ) const;

        G4double GetXshift() const {return fXshift;};
        G4double GetYshift() const {return fYshift;};
        G4double GetOrzrot() const {return fOrzrot;};
        G4double GetVerrot() const {return fVerrot;};

        //Setters
        //
        void SetXshift(const G4double& val) {fXshift=val;};
        void SetYshift(const G4double& val) {fYshift=val;};
        void SetOrzrot(const G4double& val) {fOrzrot=val;};
        void SetVerrot(const G4double& val) {fVerrot=val;};
        void SetFiberMapPath(const G4String& val) {fFiberMapPath=val;};
        G4String GetFiberMapPath() const {return fFiberMapPath;};
        void SetTowerMapPath(const G4String& val) {fTowerMapPath=val;};
        G4String GetTowerMapPath() const {return fTowerMapPath;};




    //
	//  Build contour in x-y plane of a module as 
	//  an hexcell shape
	//
	std::vector<G4TwoVector> calcmod(double radius, int nrow, int ncol); 

    private:

        //Mandatory method for Geant4
        //
        G4VPhysicalVolume* DefineVolumes();

        //Per-fiber geometry map (CSV) support
        //
        //Recorded during DefineVolumes(): one module placement per active tower,
        //and the fiber-local table (placed once into the shared module volume).
        struct ModulePlacement { G4int towerID; G4double mx, my; G4int gridCol, gridRow; };
        struct FiberLocal { char type; G4int colInMod, rowInMod, copyNumber;
                            G4double lx, ly; std::string pvName; };
        std::vector<ModulePlacement> fModulePlacements;
        std::vector<FiberLocal> fFiberLocals;
        G4double fTubeRadius{1.*mm}; //tube radius, for row/column binning tolerance
        G4double fModuleX{0.}, fModuleY{0.}; //module footprint (width,height), for tower map
        //Write the per-fiber CSV to fFiberMapPath (called once, master thread).
        void WriteFiberMap() const;
        //Write the per-tower CSV to fTowerMapPath (called once, master thread).
        void WriteTowerMap() const;

        //Members
        //
        G4bool  fCheckOverlaps; // option for checking volumes overlaps
				
				//G4VPhysicalVolume* fLeakCntPV; //PV: lekage counter
                G4VPhysicalVolume* fLeakCntlPV; //PV: lateral leakage counter
                G4VPhysicalVolume* fLeakCntdPV; //PV: longitudinal leakage counter
				G4VPhysicalVolume* fWorldPV;   //PV: wourld volume

        //Pointer to messenger for UI
        //
        DREMTubesGeoMessenger* fGeoMessenger = nullptr;

        //Parameters selectable via UI
        //
        G4double fXshift{0.}, fYshift{0.}, fVerrot{0.}, fOrzrot{0.};

        //Paths for the geometry maps; empty = disabled
        //
        G4String fFiberMapPath{""};
        G4String fTowerMapPath{""};

};

inline G4int DREMTubesDetectorConstruction::GetTowerID( const G4int& cpno ) const {
// remap as for 2021 hardware numbering from front face
//    const G4int idmap[9]={1,2,3,4,0,5,6,7,8};
// test:remap as for output of old simulation
//    const G4int idmap[9]={3,2,1,5,0,4,8,7,6};
//    return idmap[cpno-1];
    return cpno;		
}

inline G4int DREMTubesDetectorConstruction::GetSiPMID( const G4int& cpno ) const {
// kept for compatibility with old simulation. Dummy for now
    return cpno;		
}

inline G4int DREMTubesDetectorConstruction::GetSiPMTower( const G4int& town ) const {
    G4int SiPMTower=-1;
    for(int i=0;i<NoModulesSiPM;i++){
      if(town==SiPMMod[i])SiPMTower=i;
    }
    return SiPMTower;		
}

//inline const G4VPhysicalVolume* DREMTubesDetectorConstruction::GetLeakCntPV() const {
//    return fLeakCntPV;
//}
inline const G4VPhysicalVolume* DREMTubesDetectorConstruction::GetLeakCntlPV() const {
    return fLeakCntlPV;
}
inline const G4VPhysicalVolume* DREMTubesDetectorConstruction::GetLeakCntdPV() const {
    return fLeakCntdPV;
}

inline const G4VPhysicalVolume* DREMTubesDetectorConstruction::GetWorldPV() const {
    return fWorldPV;
}

#endif

//**************************************************
