//**************************************************
// \file DREMTubesEventAction.cc
// \brief: Implementation of DREMTubesEventAction 
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DREMTubesEventAction.hh"
#include "DREMTubesRunAction.hh"
#include "DREMTubesDetectorConstruction.hh"
//Includers from Geant4
//
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

//Includers from C++
//
#include <iomanip>
#include <vector>
#include <numeric>


//Define constructor
//
DREMTubesEventAction::DREMTubesEventAction()
    : G4UserEventAction(),
    EnergyScin(0.),
    EnergyCher(0.),
    NofCherDet(0),
    NofScinDet(0),
    EnergyTot(0.),
    PrimaryPDGID(0),
    PrimaryX(0),
    PrimaryY(0),
    PrimaryParticleEnergy(0.),
    //EscapedEnergy(0.),
    EscapedEnergyl(0.),
    EscapedEnergyd(0.),    
    PSEnergy(0.),
    VertexX(0.),
    VertexY(0.),
    ImpactX(0.),
    ImpactY(0.),
    VectorSignals(0.),
    VectorSignalsCher(0.),
    VecSPMT(0.),
    VecCPMT(0.),
    VecTowerE(0.) {
}

//Define de-constructor
//
DREMTubesEventAction::~DREMTubesEventAction() {}

//Define BeginOfEventAction() and EndOfEventAction() methods
//
void DREMTubesEventAction::BeginOfEventAction(const G4Event*) {  
    
    //Initialize data memebers at begin of each event
    //
    EnergyScin = 0.;
    EnergyCher = 0.;
    NofCherDet = 0;
    NofScinDet = 0;
    EnergyTot = 0.;
    PrimaryPDGID = 0;
    PrimaryX = 0;
    PrimaryY = 0;
    PrimaryParticleEnergy = 0.;
    //EscapedEnergy = 0.;
    EscapedEnergyl = 0.;
    EscapedEnergyd = 0.;    
    PSEnergy = 0.;

    VectorSignals.clear();
    VectorSignalsCher.clear();
    VecSPMT.clear();
    VecCPMT.clear();
    VecTowerE.clear();
    VecLeakCounter.clear();

    VectorSignals.assign(NoFibersTower*NoModulesSiPM, 0.);
    VectorSignalsCher.assign(NoFibersTower*NoModulesSiPM, 0.);
    VecSPMT.assign(NoModulesActive, 0.);
    VecCPMT.assign(NoModulesActive, 0.);
    VecTowerE.assign(NoModulesActive, 0.);
    VecLeakCounter.assign(4*NofLeakCounterLayers+1, 0.);


}

// Straight-line extrapolation of the primary from its vertex to the calorimeter
// front face. The front face is the plane z = -halfZ of the "Calorimeter" box,
// placed (unrotated) inside "CalorimeterBox", which carries the platform
// rotation and shift (/tbgeo/). Result in the calorimeter frame: module centres
// are at x = -128 mm*(column-2), y = 27.728 mm*(row-(NofmodulesY-1)/2) (0-based).
void DREMTubesEventAction::ComputeImpact(const G4Event* event) {
    VertexX = VertexY = ImpactX = ImpactY = 0.;
    const G4PrimaryVertex* vtx = event->GetPrimaryVertex(0);
    if (!vtx || !vtx->GetPrimary(0)) return;
    G4ThreeVector pos = vtx->GetPosition();
    G4ThreeVector dir = vtx->GetPrimary(0)->GetMomentumDirection();
    VertexX = pos.x();
    VertexY = pos.y();
    static G4ThreadLocal const G4VPhysicalVolume* boxPV = nullptr;
    static G4ThreadLocal G4double halfZ = 0.;
    if (!boxPV) {
        auto store = G4PhysicalVolumeStore::GetInstance();
        boxPV = store->GetVolume("CalorimeterBox", false);
        auto caloPV = store->GetVolume("Calorimeter", false);
        if (caloPV) halfZ = static_cast<const G4Box*>(caloPV->GetLogicalVolume()->GetSolid())->GetZHalfLength();
        if (!boxPV || !caloPV) { G4cerr << "ComputeImpact: calorimeter volumes not found" << G4endl; return; }
    }
    G4RotationMatrix rinv = boxPV->GetObjectRotationValue().inverse();
    G4ThreeVector p = rinv * (pos - boxPV->GetObjectTranslation());
    G4ThreeVector d = rinv * dir;
    if (d.z() == 0.) return;
    G4double t = (-halfZ - p.z()) / d.z();
    ImpactX = p.x() + t * d.x();
    ImpactY = p.y() + t * d.y();
}

void DREMTubesEventAction::EndOfEventAction(const G4Event* event) {
 
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    ComputeImpact(event);

    //Add all p.e. in Scin and Cher fibers before calibration
    //
    NofScinDet = std::accumulate(VecSPMT.begin(),VecSPMT.end(),0);
    NofCherDet = std::accumulate(VecCPMT.begin(),VecCPMT.end(),0);

    G4int NofPheSciSiPM = std::accumulate(VectorSignals.begin(),VectorSignals.end(),0.);
    G4int NofPheCerSiPM = std::accumulate(VectorSignalsCher.begin(),VectorSignalsCher.end(),0.);

    // NofScinDet and NofCherDet take sum of photoelectrons deposited in PMT towers
    // also those in the SiPM modules, as if they were readout with PMTs
    // to add: sum separately photoelectrons in PMT and SiPM towers

    //Fill ntuple event by event
    //entries with vectors are automatically filled
    //
    analysisManager->FillNtupleDColumn(0, EnergyScin);
    analysisManager->FillNtupleDColumn(1, EnergyCher);
    analysisManager->FillNtupleDColumn(2, NofCherDet);
    analysisManager->FillNtupleDColumn(3, NofScinDet);
    analysisManager->FillNtupleDColumn(4, EnergyTot);
    analysisManager->FillNtupleDColumn(5, PrimaryParticleEnergy);
    analysisManager->FillNtupleIColumn(6, PrimaryPDGID);
    analysisManager->FillNtupleDColumn(7, EscapedEnergyl);
    analysisManager->FillNtupleDColumn(8, EscapedEnergyd);
    analysisManager->FillNtupleDColumn(9, PSEnergy);
    analysisManager->FillNtupleDColumn(10, PrimaryX);
    analysisManager->FillNtupleDColumn(11,PrimaryY);
    analysisManager->FillNtupleDColumn(12,NofPheSciSiPM);
    analysisManager->FillNtupleDColumn(13,NofPheCerSiPM);    
    analysisManager->FillNtupleDColumn(14,VertexX);
    analysisManager->FillNtupleDColumn(15,VertexY);
    analysisManager->FillNtupleDColumn(16,ImpactX);
    analysisManager->FillNtupleDColumn(17,ImpactY);
    analysisManager->AddNtupleRow();
    //Vector entries in ntuple are automatically filled

}

//**************************************************
