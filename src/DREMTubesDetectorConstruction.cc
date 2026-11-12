//**************************************************
// \file DREMTubesDetectorConstruction.cc
// \brief: Implementation of 
//         DREMTubesDetectorConstruction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DREMTubesDetectorConstruction.hh"
#include "DREMTubesGeoMessenger.hh"

//Includers from Geant4
//
#include <random>
#include <iostream>
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Sphere.hh"
#include "G4Colour.hh"
#include "G4TwoVector.hh"

//Messenger constructor
//
DREMTubesGeoMessenger::DREMTubesGeoMessenger(DREMTubesDetectorConstruction* DetConstruction)
    : fDetConstruction(DetConstruction){
    
    //The Messenger directory and commands must be initialized
    //within constructor
    //
    fMsgrDirectory = new G4UIdirectory("/tbgeo/");
    fMsgrDirectory->SetGuidance("Set movable parameters in test-beam geometry.");

    fXshiftcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/xshift", this);
    fXshiftcmd->SetParameterName("xshift",/*omittable=*/false,/*currentAsDefault=*/true);
    fXshiftcmd->SetGuidance("Shift test-beam platform x direction (default unit mm)");
    fXshiftcmd->SetDefaultUnit("mm");
    fXshiftcmd->SetDefaultValue(0.);
    fYshiftcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/yshift", this);
    fYshiftcmd->SetParameterName("yshift",/*omittable=*/false,/*currentAsDefault=*/true);
    fYshiftcmd->SetGuidance("Shift test-beam platform y direction (default unit mm)");
    fYshiftcmd->SetDefaultUnit("mm");
    fYshiftcmd->SetDefaultValue(0.);
    fOrzrotcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/horizrot", this);
    fOrzrotcmd->SetParameterName("horizrot",/*omittable=*/false,/*currentAsDefault=*/true);
    fOrzrotcmd->SetGuidance("Rotate platform (default deg)");
    fOrzrotcmd->SetDefaultUnit("deg");
    fOrzrotcmd->SetDefaultValue(0.);
    fVerrotcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/vertrot", this);
    fVerrotcmd->SetParameterName("vertrot",/*omittable=*/false,/*currentAsDefault=*/true);
    fVerrotcmd->SetGuidance("Lift up calorimeter from back side (default deg)");
    fVerrotcmd->SetDefaultUnit("deg");
    fVerrotcmd->SetDefaultValue(0.);
}


//Messenger destructor
//
DREMTubesGeoMessenger::~DREMTubesGeoMessenger(){

    //The Messenger fields should be deleted
    //in destructor
    delete fMsgrDirectory;
    delete fXshiftcmd;
    delete fYshiftcmd;
    delete fOrzrotcmd;
    delete fVerrotcmd;
}

//Messenger SetNewValue virtual method from base class
//
void DREMTubesGeoMessenger::SetNewValue(G4UIcommand* command, G4String newValue){

    if(command == fXshiftcmd){
        fDetConstruction->SetXshift(fXshiftcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: x-shifted test-beam setup by "<<fDetConstruction->GetXshift()<<" mm"<<G4endl;
    }
    else if(command == fYshiftcmd){
        fDetConstruction->SetYshift(fYshiftcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: y-shifted test-beam setup by "<<fDetConstruction->GetYshift()<<" mm"<<G4endl;
    }
    else if(command == fOrzrotcmd){
        fDetConstruction->SetOrzrot(fOrzrotcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: orz-rotated test-beam setup by "<<fDetConstruction->GetOrzrot()<<" rad"<<G4endl;
    }
    else if(command == fVerrotcmd){
        fDetConstruction->SetVerrot(fVerrotcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: ver-rotated test-beam setup by "<<fDetConstruction->GetVerrot()<<" rad"<<G4endl;
    }
}


//
//  sqrt3 constants used in code
//  reciprocal of sqrt3 given as number that divided 
//  by sqrt 3 gives 3 
//
const G4double sq3=1.733;
const G4double sq3m1=sq3/3.;

//Constructor
//
DREMTubesDetectorConstruction::DREMTubesDetectorConstruction()
    : G4VUserDetectorConstruction(),
    fCheckOverlaps(false),
		//fLeakCntPV(nullptr),
    fLeakCntlPV(nullptr),
    fLeakCntdPV(nullptr),
    fWorldPV(nullptr){

    fGeoMessenger = new DREMTubesGeoMessenger(this);

}

//De-constructor
//
DREMTubesDetectorConstruction::~DREMTubesDetectorConstruction() {    delete fGeoMessenger;}

//Define Construct() method
G4VPhysicalVolume* DREMTubesDetectorConstruction::Construct() {
  
    // Define volumes
    return DefineVolumes();
}

G4VPhysicalVolume* DREMTubesDetectorConstruction::DefineVolumes() {

    //--------------------------------------------------
    //Define Elements, Mixtures and Materials
    //--------------------------------------------------

    auto nistManager = G4NistManager::Instance();

    //Elements
    //
    G4String name, symbol;    
    G4double a, z;            // a=mass of a mole, z=mean number of protons;  
  
    a = 1.01*g/mole;
    G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a); //Hidrogen

    a = 12.01*g/mole;
    G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a); //Carbon

    a = 16.00*g/mole;
    G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a); //Oxygen

    a = 28.09*g/mole;
    G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a); //Silicon
  
    a = 18.9984*g/mole;
    G4Element* elF  = new G4Element("Fluorine",symbol="F" , z= 9., a); //Fluorine

    //a = 63.546*g/mole;
    //G4Element* elCu = new G4Element("Copper", symbol="Cu", z=29., a); //Copper
    auto elCu = nistManager->FindOrBuildElement(29, true);

    //a = 65.38*g/mole;
    //G4Element* elZn = new G4Element("Zinc", symbol="Zn", z=30., a); //Zinc
    auto elZn = nistManager->FindOrBuildElement(30, true);
    auto elFe = nistManager->FindOrBuildElement(26, true);
    auto elCr = nistManager->FindOrBuildElement(24, true);
    auto elNi = nistManager->FindOrBuildElement(28, true);

    //Materials 
    //

    // Polystyrene from elements (C5H5)
    G4Material* Polystyrene = new G4Material("Polystyrene", 1.05*g/cm3, 2);
    Polystyrene->AddElement(elC, 8);
    Polystyrene->AddElement(elH, 8); 

    // PMMA material from elements (C502H8)
    // 
    auto PMMA = new G4Material("PMMA", 1.19*g/cm3, 3); 
    PMMA->AddElement(elC, 5);
    PMMA->AddElement(elO, 2);
    PMMA->AddElement(elH, 8); 
    
    // Fluorinated Polymer material from elements (C2F2)
    // material for the cladding of the Cherenkov fibers
    auto fluorinatedPolymer = new G4Material("Fluorinated_Polymer", 1.43*g/cm3, 2);
    fluorinatedPolymer->AddElement(elC,2);
    fluorinatedPolymer->AddElement(elF,2);

    // Glass material from elements (SiO2)
    //
    auto Glass = new G4Material("Glass", 2.4*g/cm3, 2);
    Glass -> AddElement(elSi, 1);
    Glass -> AddElement(elO, 2); 

    // Mixtures
    //

    // CuZn37 (Brass)
    //
    const double BrassDensity = 8.44*g/cm3;
    auto CuZn37 = new G4Material(name="Brass", BrassDensity, 2);
    CuZn37->AddElement(elCu, 0.7);
    CuZn37->AddElement(elZn, 0.3);

    // 304 stainless steel 
    auto Steel304 = new G4Material("Steel304", 7.93*g/cm3, 3);
    Steel304->AddElement(elFe, 0.68);
    Steel304->AddElement(elCr, 0.20);
    Steel304->AddElement(elNi, 0.12);

    // Assign material to the calorimeter volumes
    //
    G4Material* defaultMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
    //G4Material* absorberMaterial = nistManager->FindOrBuildMaterial("G4_Cu");
    G4Material* SiMaterial = nistManager->FindOrBuildMaterial("G4_Si");
    G4Material* LeadMaterial = nistManager->FindOrBuildMaterial("G4_Pb");
    G4Material* PSScinMaterial = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");
    //G4Material* absorberMaterial = G4Material::GetMaterial("Brass");
    G4Material* absorberMaterial = G4Material::GetMaterial("Steel304");
    G4Material* ScinMaterial = G4Material::GetMaterial("Polystyrene");
    G4Material* CherMaterial = G4Material::GetMaterial("PMMA");
    G4Material* GlassMaterial = G4Material::GetMaterial("Glass");
    G4Material* CladCherMaterial = G4Material::GetMaterial("Fluorinated_Polymer");
    G4Material* LeakCounterScinMaterial = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");

    auto invisibleAttr = new G4VisAttributes();
    invisibleAttr->SetVisibility(false);

    //--------------------------------------------------
    //Define Optical Properties
    //--------------------------------------------------

    // Use Energy(eV)=1.24/waevelenght(um)
    // 2.034eV is 610nm RED 
    // 2.75eV is 450nm BLUE (peak of scintillating fibers)
    // 3.09eV is 400nm VIOLET (end of visible)
    //4.1eV is 300nm UV (cherenkov peak is 310-350nm)
    //
    const G4int ENTRIES = 32;
    G4double photonEnergy[ENTRIES] =                    
        { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,   
          2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     
          2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
          2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
          2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV, 
          3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV, 
          3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
          3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV }; 
    G4double rindexScin[ENTRIES] =
        { 1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59 };
    /*G4double absorptionScin[ENTRIES] =
        { 400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm };*/

    G4MaterialPropertiesTable *MPTScin = new G4MaterialPropertiesTable();
    MPTScin -> AddProperty("RINDEX", photonEnergy, rindexScin, true);
    MPTScin->AddConstProperty("SCINTILLATIONYIELD", 1.0, true);   // to geant4.11
    MPTScin->AddConstProperty("EXCITATIONRATIO", 0.0, true);

    /*MPTScin -> AddProperty("ABSLENGTH",
         photonEnergy, absorptionScin, ENTRIES)->SetSpline(true);*/

    G4double rindexCher[ENTRIES] =
        { 1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49 };
    /*G4double absorptionCher[ENTRIES] = 
        { 890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm };*/

    G4MaterialPropertiesTable *MPTCher = new G4MaterialPropertiesTable();
    MPTCher -> AddProperty("RINDEX", photonEnergy, rindexCher, ENTRIES, true);
    /*MPTCher -> AddProperty("ABSLENGTH", 
            photonEnergy, absorptionCher, ENTRIES)->SetSpline(true);*/
    CherMaterial -> SetMaterialPropertiesTable(MPTCher);

    G4double rindexCherclad[ENTRIES] =
        { 1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42 };

    G4MaterialPropertiesTable *MPTCherclad = new G4MaterialPropertiesTable();
    MPTCherclad -> AddProperty("RINDEX", photonEnergy, rindexCherclad, ENTRIES, true);
    CladCherMaterial -> SetMaterialPropertiesTable(MPTCherclad);

    G4double rindexglass[ENTRIES] =
        { 1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51 };

    G4MaterialPropertiesTable *MPTglass = new G4MaterialPropertiesTable();
    MPTglass -> AddProperty("RINDEX", photonEnergy, rindexglass, ENTRIES, true);
    GlassMaterial -> SetMaterialPropertiesTable(MPTglass);

    G4double rindexSi[ENTRIES] =
        { 3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42 };

    G4double absorptionSi[ENTRIES] = 
        { 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm };

    G4MaterialPropertiesTable *MPTSi = new G4MaterialPropertiesTable();
    MPTSi -> AddProperty("RINDEX", photonEnergy, rindexSi, ENTRIES, true);
    MPTSi -> AddProperty("ABSLENGHT", photonEnergy, absorptionSi, ENTRIES, true);
    SiMaterial -> SetMaterialPropertiesTable(MPTSi); 
  
    // Scintillating proprieties of the scintillating fiber material
    // Birks constant of the polystyrene
    //
    G4double Scin_FAST[ENTRIES] = // Emission spectrum for the fast component 
        { 0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.1,
          0.2, 0.4, 0.6, 0.8,
          1., 0.8, 0.6, 0.1,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0. };
    G4double Scin_SLOW[ENTRIES] = // Emission spectrum for the slow component
        { 0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0. };

    ScinMaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

    MPTScin -> AddProperty("FASTCOMPONENT", photonEnergy, Scin_FAST, ENTRIES, true);
    MPTScin -> AddProperty("SLOWCOMPONENT", photonEnergy, Scin_SLOW, ENTRIES, true);
    MPTScin -> AddConstProperty("SCINTILLATIONYIELD", 10000./MeV, true); 
    // Typical is 10000./MeV (this is what makes full simulations long as hell)
    MPTScin->AddConstProperty("EXCITATIONRATIO", 0.0, true);
    MPTScin -> AddConstProperty("RESOLUTIONSCALE", 1.0, true); 
    // Broad the fluctuation of photons produced
    MPTScin -> AddConstProperty("FASTTIMECONSTANT", 2.8*ns, true);
    MPTScin -> AddConstProperty("SLOWTIMECONSTANT", 10.*ns, true);
    MPTScin -> AddConstProperty("YIELDRATIO", 1.0, true); 
    // I don't want a slow component, if you want it must change
    ScinMaterial -> SetMaterialPropertiesTable(MPTScin);

    //--------------------------------------------------
    //Define Volumes
    //--------------------------------------------------
    
//  tubes are packed in Y, so Y distance between two tubes is sqrt(3)/2*2.*tuberadius
    double tolerance = 0.0*mm;
    G4double tuberadius = 1.0*mm;
    G4double dtubeY=sq3*tuberadius;
    G4double dtubeX=2.*tuberadius;
    G4double moduleX = (2*NofFiberscolumn+1)*tuberadius; 
    G4double moduleY = (NofFibersrow-1.)*dtubeY+4.*tuberadius*sq3m1;
//  side of hexagon in which tube is inscribed
//    G4double side=tuberadius*2*sq3m1;

    // Geometry parameters of the world, world is a G4Box
    //
    //G4double worldX = 200 * moduleX;
    //G4double worldY = 200 * moduleY;
    //G4double worldZ = 60 * moduleZ;

    // Geometry parameters of the fiber
    //
    //G4double fiberradius = 0.5*mm;
    G4double fiberZ = moduleZ;
    
    // Geometry parameters of the core
    //
    G4double coreradius = 0.485*mm;
    G4double coreZ = moduleZ;

    // Geometry parameters of the cladding
    //
    G4double claddingradiusmin = 0.485*mm;
    G4double claddingradiusmax = 0.50*mm;
    G4double claddingZ = moduleZ;
    
    // Geometry parameters of the tube
    //
    //G4double tubeZ = fiberZ;

    // Geometry parameters of the SiPM
    //
    G4double SiPMX = 1.*mm;
    G4double SiPMY = SiPMX;
    G4double SiPMZ = 0.36*mm;
    G4double SiPMZh = SiPMZ/2.;

    // Geometry parameters of the SiPM, active silicon layer
    //
    G4double SiX = 1.*mm;
    G4double SiY = SiX;
    G4double SiZ = 0.05*mm;

    // Geometry parameters of the module equipped with SiPM
    //
    G4double moduleequippedZ = moduleZ + SiPMZ;
    G4double moduleequippedX = moduleX; 
    G4double moduleequippedY = moduleY;

    //Preshower dimensions
    //
    G4double PSX = 9.2*cm;
    G4double PSY = 9.2*cm;
    G4double PSZ = 1.*cm;


    // Leakage counter dimensions
    G4double leakBoxX = 50.*cm; 
    G4double leakBoxY = 12.*cm;
    G4double leakBoxZ = 50.*cm;    
    G4double tailCatcherDist = 50.*cm;


    // Define calorimeter size, and then define world on it
    // to ensure world is large enough
    // Calorimeter (matrix of modules equipped)
    // 
    G4double caloX=0.;
    G4double caloY=0.;
    G4double caloZ=0.;
    if(irot) {
      caloX=moduleequippedY*NofmodulesX/2;
      caloY=moduleequippedX*NofmodulesY/2;
      caloZ=moduleequippedZ/2;      
    }
    else {
      caloX=moduleequippedX*NofmodulesX/2;
      caloY=moduleequippedY*NofmodulesY/2;
      caloZ=moduleequippedZ/2;      
    }
    G4double worldX = 50. * caloX;
    G4double worldY = 50. * caloY;
    G4double worldZ = 50. * caloZ;
    std::cout << " caloX " << caloX << " caloY " << caloY << std::endl;
   

    // Building geometries
    //
    // World
    //
    G4VSolid* worldS  = new G4Box("World", worldX/2, worldY/2, worldZ/2); 
                         
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,          
                                                   defaultMaterial, 
                                                   "World");       
  
    //worldLV->SetVisAttributes(G4VisAttributes::Invisible);
    worldLV->SetVisAttributes(invisibleAttr);


    fWorldPV = new G4PVPlacement( 0,                // no rotation
                                  G4ThreeVector(),  // at (0,0,0)
                                  worldLV,          // its logical
                                  "World",          // its name
                                  0,                // its mother
                                  false,            // no boolean oper 
                                  0,                // copy number
                                  fCheckOverlaps);  // check overlaps 

    //Preshower
    //
    if(PreShowerIn)
    {
      auto PSSolid = new G4Box("Preshower", PSX/2., PSY/2., PSZ/2.);

      auto PSLV = new G4LogicalVolume(PSSolid, defaultMaterial, "Preshower");

      new G4PVPlacement( 0, 
            //G4ThreeVector(0.,0.,-335.*cm),
            //G4ThreeVector(0.,0.,-moduleZ/2 - 15.5*cm),
            // This should be corrected with TB setup; note that PS position should be outside calo box, otherwise data is not saved           
            //G4ThreeVector(0.,0.,-caloZ - 60*cm), 
            //G4ThreeVector(0.,0.,-caloZ - 150.*cm), 
            G4ThreeVector(0.,0.,-caloZ - 70.*cm), 

            PSLV,
            "Preshower",
            worldLV,
            false,
            0,
            fCheckOverlaps);	 

      auto PSLeadSolid = new G4Box("Preshower_pb", PSX/2., PSY/2., PSZ/4.);

      auto PSLeadLV = new G4LogicalVolume(PSLeadSolid, LeadMaterial, "Preshower_pb");

      new G4PVPlacement( 0, 
            G4ThreeVector(0.,0.,-PSZ/4.),
            PSLeadLV,
            "Preshower_pb",
            PSLV,
            false,
            0,
            fCheckOverlaps);	 

      G4VisAttributes* PbVisAtt = new G4VisAttributes( G4Colour::Grey() );
      PbVisAtt->SetVisibility(true);
      PbVisAtt->SetForceSolid(true);
      PSLeadLV->SetVisAttributes( PbVisAtt );

      auto PSScinSolid = new G4Box("Preshower_scin", PSX/2., PSY/2., PSZ/4.);

      auto PSScinLV = new G4LogicalVolume(PSScinSolid, PSScinMaterial, "Preshower_scin");

      new G4PVPlacement( 0, 
            G4ThreeVector(0.,0.,PSZ/4.),
                        PSScinLV,
                  "Preshower_scin",
                        PSLV,
                        false,	
                        0,
                        fCheckOverlaps);	 

      G4VisAttributes* PSScinVisAtt = new G4VisAttributes( G4Colour::Cyan() );
      PSScinVisAtt->SetVisibility(true);
      PSScinLV->SetVisAttributes( PSScinVisAtt );
      // End Preshower
    }





    /********** Original leakage volume (sphere)  ***************/
    /*
    //Absorber to calculate leakage
    //
    G4VSolid* leakageabsorber = new G4Sphere("leakageabsorber",                        
        17.*m, 17.1*m, 0.*deg, 360.*deg, 0.*deg, 180.*deg); 
    
    G4LogicalVolume* leakageabsorberLV = new G4LogicalVolume(leakageabsorber,
                                                             defaultMaterial,  
                                                             "leakageabsorber");        
    
    leakageabsorberLV->SetVisAttributes(G4VisAttributes::Invisible);   

    fLeakCntPV = new G4PVPlacement( 0, G4ThreeVector(),
				    leakageabsorberLV,         
                                    "leakageabsorber",
                                    worldLV,               
                                    false,          
                                    0,               
                                    fCheckOverlaps);
    */

                                    
   // Module equipped (with SiPM)
   //
   // Basic module structure: extrusion of an hexcell shape
    G4TwoVector offA(0,0), offB(0,0);
    G4double scaleA = 1, scaleB = 1;
    auto polygon=calcmod(tuberadius, NofFibersrow, NofFiberscolumn);
    G4VSolid* moduleequippedS = new G4ExtrudedSolid("moduleequipped", 
		                                     polygon, 
						     moduleequippedZ/2, 
						     offA, scaleA, offB, scaleB);
    G4LogicalVolume* moduleequippedLV = new G4LogicalVolume(moduleequippedS,
                                                            defaultMaterial,
                                                            "moduleequipped"); 


    // Define Calorimeter Box to put physical leakage counters in
    // distances from calorimeter to box surface ~ 10cm
    G4double distCaloBox = 5.*cm;
    G4double caloBoxX = caloX+leakBoxX/2+distCaloBox; G4double caloBoxY = caloY+leakBoxZ/2+distCaloBox; G4double caloBoxZ = caloZ + distCaloBox;
    //G4VSolid* CalorimeterBoxSolid = new G4Box("CalorimeterBoxSolid", caloX+leakBoxX/2+distCaloBox, caloY+leakBoxZ/2+distCaloBox, caloZ);
    G4VSolid* CalorimeterBoxSolid = new G4Box("CalorimeterBoxSolid", caloBoxX, caloBoxY, caloBoxZ);

    G4LogicalVolume* CalorimeterBoxLV = new G4LogicalVolume( CalorimeterBoxSolid,
                                                          defaultMaterial,
                                                          "CalorimeterBoxLV");



    G4VSolid* CalorimeterS = new G4Box("CalorimeterS",caloX,caloY,caloZ);

    G4LogicalVolume* CalorimeterLV = new G4LogicalVolume( CalorimeterS,
                                                          defaultMaterial,
                                                          "CalorimeterLV");

    // Modules equipped placement
    //
    G4int copynumbermodule = 0;
    G4double m_x, m_y;
    G4ThreeVector vec_m;
    for(int column=0; column<NofmodulesX; column++){ 
        for(int row=0; row<NofmodulesY; row++){
	    G4int ii=column+row*NofmodulesX;
            if(irot) {
              m_x = -dtubeY*NofFibersrow*column+ dtubeY*NofFibersrow*(NofmodulesX-1)/2;
              m_y = row*NofFiberscolumn*dtubeX-NofFiberscolumn*dtubeX*(NofmodulesY-1)/2;
            }
            else {
              m_x = -dtubeX*NofFiberscolumn*column+ dtubeX*NofFiberscolumn*(NofmodulesX-1)/2;
              m_y = row*NofFibersrow*dtubeY-NofFibersrow*dtubeY*(NofmodulesY-1)/2;
            }	     	    
            if(modflag[ii]>=0) {        
              copynumbermodule = modflag[ii];
//              copynumbermodule = (1+column)+(row*NofmodulesX);
//	      std::cout << " column " << column << " row " << row << " cpnm " << copynumbermodule << std::endl;
// setup for 90 deg rotation of modules
              if(irot){
                G4RotationMatrix rotmod  = G4RotationMatrix();
                rotmod.rotateZ(90.*degree);
                G4ThreeVector posm;
                posm.setX(m_x);
                posm.setY(m_y);
                posm.setZ(0.);
                G4Transform3D transfm = G4Transform3D(rotmod,posm);
                new G4PVPlacement(transfm,
                                  moduleequippedLV,     
                                  "moduleequipped",                        
                                  CalorimeterLV,                      
                                  false,                          
                                  copynumbermodule); 
	      }
	      else {
// simple positioning
                vec_m.setX(m_x);
                vec_m.setY(m_y);
                vec_m.setZ(0.);
                new G4PVPlacement(0,
                                  vec_m,              
                                  moduleequippedLV,     
                                  "moduleequipped",                        
                                  CalorimeterLV,                      
                                  false,                          
                                  copynumbermodule); 
              } // irot
	    } //ifmod
        };
    }; 
 
    // Calorimeter placement (with rotation wrt beam axis)
    //
    // Negative sign to make angles same as TB platform
    G4double xrot = - fVerrot;
    G4double yrot = - fOrzrot;
    G4RotationMatrix rotX = G4RotationMatrix();
    G4RotationMatrix rotY = G4RotationMatrix();
    rotX.rotateX(xrot);
    rotY.rotateY(yrot);
    G4RotationMatrix rotm = rotY * rotX;

    G4ThreeVector position;

    position.setX(fXshift);
    position.setY(fYshift);
    position.setZ(0.);

    G4Transform3D transform = G4Transform3D(rotm,position); 


    //
    //  Build closed tube for detailed leakage study
    //
    if(TruthLeakageIn)
    {
      G4double leakradint=sqrt( (caloBoxX)*(caloBoxX)+(caloBoxY)*(caloBoxY));	// Added *1.2 wrt Giacomo's

      G4double leakradout=leakradint+20*cm;
      G4double tube_dPhi = 2.* M_PI * rad;
      G4double disc_th = 20.*cm;
      G4VSolid* leakageabsorberl = new G4Tubs("leakageabsorberl",
          leakradint, leakradout, caloBoxZ + tailCatcherDist + leakBoxZ, 0., tube_dPhi );

      G4LogicalVolume* leakageabsorberlLV = new G4LogicalVolume(leakageabsorberl,
                                                              defaultMaterial,  
                                                              "leakageabsorberl");        
      G4VisAttributes* LkVisAttl = new G4VisAttributes(G4Colour(0.0,0.8,0.0)); //green
      LkVisAttl->SetVisibility(true);
      //LkVisAttl->SetForceWireframe(true);
      LkVisAttl->SetForceSolid(true);
      leakageabsorberlLV->SetVisAttributes(LkVisAttl);
      //leakageabsorberlLV->SetVisAttributes(invisibleAttr); 
      //leakageabsorberlLV->SetVisAttributes(G4VisAttributes::Invisible);

      fLeakCntlPV = new G4PVPlacement( transform,
              leakageabsorberlLV,         
                                      "leakageabsorberl",
                                      worldLV,               
                                      false,          
                                      0,               
                                      fCheckOverlaps);

      G4VSolid* leakageabsorberd = new G4Tubs("leakageabsorberd",
          0, leakradout, disc_th/2, 0., tube_dPhi );

      G4LogicalVolume* leakageabsorberdLV = new G4LogicalVolume(leakageabsorberd,
                                                              defaultMaterial,  
                                                              "leakageabsorberd");        
      G4VisAttributes* LkVisAttd = new G4VisAttributes(G4Colour(0.0,0.4,0.0)); //green
      LkVisAttd->SetVisibility(true);
      //LkVisAttd->SetForceWireframe(true);
      LkVisAttd->SetForceSolid(true);
      leakageabsorberdLV->SetVisAttributes(LkVisAttd);
      
      
      //leakageabsorberdLV->SetVisAttributes(invisibleAttr);
      G4ThreeVector positiond;
      positiond.setX(caloBoxZ*sin(xrot));
      positiond.setY(-caloBoxZ*sin(yrot));

      positiond.setZ(caloBoxZ+disc_th/2 + tailCatcherDist + 2*leakBoxY + 30.*cm);
      G4Transform3D transformd = G4Transform3D(rotm,positiond); 
      fLeakCntdPV = new G4PVPlacement( transformd,
                        leakageabsorberdLV,
                        "leakageabsorberd",
                        worldLV,
                        false,
                        0,
                        fCheckOverlaps);
    } // end truth leak box





    // rotate box containing calorimeter according to input angles and shifts
    /*G4VPhysicalVolume* CalorimeterPV =*/ new G4PVPlacement(transform,
                                                         CalorimeterBoxLV,
                                                         "CalorimeterBox",
                                                         worldLV,
                                                         false,
                                                         0,
                                                         fCheckOverlaps);


    // Calorimter rotation inside box (just for a correct placing)
    G4RotationMatrix RotationInsideBox  = G4RotationMatrix();
    G4ThreeVector PositionInsideBox;
    PositionInsideBox.setX(0.); PositionInsideBox.setY(0.); PositionInsideBox.setZ(0.);    
    G4Transform3D TransformInsideBox = G4Transform3D(RotationInsideBox, PositionInsideBox); 


    /*G4VPhysicalVolume* CalorimeterPV =*/ new G4PVPlacement(TransformInsideBox,
                                                         CalorimeterLV,
                                                         "Calorimeter",
                                                         CalorimeterBoxLV,
                                                         false,
                                                         0,
                                                         fCheckOverlaps);




    // Module (same shape as moduleequipped)
    //
    G4VSolid* moduleS = new G4ExtrudedSolid("module", 
		                            polygon, 
					    moduleZ/2, 
					    offA, scaleA, offB, scaleB);
                         
    G4LogicalVolume* moduleLV = new G4LogicalVolume(moduleS,
                                                    defaultMaterial,
                                                    "module");

    G4ThreeVector pos_module;
    pos_module.setX(0.);
    pos_module.setY(0.);
    pos_module.setZ(-0.18);
                              
    /*G4VPhysicalVolume* modulePV =*/ new G4PVPlacement(0,
                                                    pos_module,
                                                    moduleLV,
                                                    "module",
                                                     moduleequippedLV,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);

    // Optical Surface properties between the glass and the Si of the SiPM
    G4OpticalSurface* OpSurfaceGlassSi = new G4OpticalSurface("OpSurfaceGlassSi");
    OpSurfaceGlassSi -> SetType(dielectric_metal);
    OpSurfaceGlassSi -> SetModel(glisur);
    OpSurfaceGlassSi -> SetFinish(polished);
    G4double efficiencyOpSurfaceGlassSi[ENTRIES] =     //100% detection efficiency 
                                    { 1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1 };

    /*G4double efficiencyOpSurfaceGlassSi[ENTRIES] =     //0% detection efficiency 
                                    { 0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0 };*/

    G4double reflectivityOpSurfaceGlassSi[ENTRIES] =  // 0% reflection
                                    { 0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0. };

    G4MaterialPropertiesTable* MPTOpSurfaceGlassSi = new G4MaterialPropertiesTable();
    MPTOpSurfaceGlassSi -> AddProperty("EFFICIENCY", photonEnergy, efficiencyOpSurfaceGlassSi, ENTRIES, true);
    MPTOpSurfaceGlassSi -> AddProperty("REFLECTIVITY", photonEnergy, reflectivityOpSurfaceGlassSi, ENTRIES, true);
    OpSurfaceGlassSi -> SetMaterialPropertiesTable(MPTOpSurfaceGlassSi);

    // SiPM
    //
    G4VSolid* SiPMS = new G4Box("SiPM", SiPMX/2, SiPMY/2, SiPMZ/2);
                         
    G4LogicalVolume* SiPMLV = new G4LogicalVolume(SiPMS, GlassMaterial,"SiPM");

//    SiPMLV->SetVisAttributes(G4VisAttributes::Invisible);
    // Here I build the Si of the SiPM
    // 
    G4VSolid* SiS = new G4Box("Si", SiX/2, SiY/2, SiZ/2);
                         
    G4LogicalVolume* SiLV = new G4LogicalVolume( SiS, SiMaterial, "Si");

    // Si placement inside SiPM
    //
    G4ThreeVector vec_Si;
    vec_Si.setX(0.);
    vec_Si.setY(0.);
    vec_Si.setZ(SiPMZ/2-SiZ/2); // Si at the end of SiPM
                             
    /*G4VPhysicalVolume* SiPV =*/ new G4PVPlacement(0,
                                                vec_Si,  
                                                SiLV,
                                                "Si",
                                                SiPMLV,
                                                false,
                                                0,
                                                fCheckOverlaps);
    G4VisAttributes* SiVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0)); //green
    SiVisAtt->SetVisibility(true);
    SiVisAtt->SetForceWireframe(true);
    SiVisAtt->SetForceSolid(true);
    SiLV->SetVisAttributes(SiVisAtt);
//    SiLV->SetVisAttributes(G4VisAttributes::Invisible);

    // Logical Skin Surface placement around the silicon of the SiPM
    //
    /*G4LogicalSkinSurface* OpsurfaceSi =*/ new G4LogicalSkinSurface("OpsurfaceSi", 
        SiLV, OpSurfaceGlassSi);

    // Optical Surface properties between the scintillating fibers
    // and the default material
    // I'm trying to define an optical surface completly blacked 
    // as if we absorb the light at one end of fibers
    //
    G4OpticalSurface* OpSurfacedefault = new G4OpticalSurface("OpSurfacedefault");
    OpSurfacedefault -> SetType(dielectric_dielectric);
    OpSurfacedefault -> SetModel(unified);
    OpSurfacedefault -> SetFinish(polishedbackpainted); 
    // Painted from inside the fibers, light is absorbed

    // Tubes with scintillating fibers and SiPM next to them
    //
    // Attention: I place an optical surface painted (blacked) from the moduleequippedPV 
    // to the SiPMPV, in so doing I completly avoid any cross talk between SiPMs
    //
    //G4VPhysicalVolume* physi_S_fiber[NofFiberscolumn][NofFibersrow];
    //G4VPhysicalVolume* physi_SiPM[NofFiberscolumn][NofFibersrow];  
    //G4LogicalBorderSurface* logic_OpSurface_defaultAir[NofFiberscolumn][NofFibersrow];
		
    G4int copynumber = 0;

    for(int column=0; column<NofFiberscolumn; column++){
        
        std::stringstream S_fiber_column;
        S_fiber_column.str("");
        S_fiber_column << column;

        for(int row=0; row<NofFibersrow; row++){
            
            std::stringstream S_fiber_row;
            S_fiber_row.str("");
            S_fiber_row << row;
            std::string S_name;
            std::string SiPM_name;
            S_name = "S_column_" + S_fiber_column.str() + "_row_" + S_fiber_row.str(); 
            SiPM_name = "S_SiPM"; 
            //SiPM_name = "SiPMS_column" + S_fiber_column.str() + "_row_" + S_fiber_row.str();

            G4double S_x, S_y;
            G4ThreeVector vec_S_fiber;
            G4ThreeVector vec_SiPM;

            if(row%2==0){
                S_x = +moduleX/2 - tuberadius - (tuberadius*2+2*tolerance)*column;
                S_y = -moduleY/2 + tuberadius + (sq3*tuberadius+2*tolerance*mm)*(row)+tuberadius*(2.*sq3m1-1.);

                vec_S_fiber.setX(S_x);
                vec_S_fiber.setY(S_y);
                vec_S_fiber.setZ(0.);

                vec_SiPM.setX(S_x);
                vec_SiPM.setY(S_y);
                vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);
            
                copynumber = ((NofFibersrow/2)*column+row/2);
                auto logic_S_fiber = constructscinfiber(tolerance,
                                                        tuberadius,
                                                        fiberZ,
                                                        absorberMaterial,
                                                        coreradius,
                                                        coreZ,
                                                        ScinMaterial,
                                                        claddingradiusmin,
                                                        claddingradiusmax,
                                                        claddingZ,
                                                        CherMaterial);
                // Tubes with scintillating fiber placement
                //
                /*physi_S_fiber[column][row] =*/ new G4PVPlacement(0,
                                                               vec_S_fiber,
                                                               logic_S_fiber,
                                                               S_name,
                                                               moduleLV,
                                                               false,
                                                               copynumber); 

                // SiPM placement
                //
                /*physi_SiPM[column][row] =*/ new G4PVPlacement(0,
                                                            vec_SiPM,
                                                            SiPMLV,
                                                            SiPM_name,
                                                            moduleequippedLV,
                                                            false,
                                                            copynumber); //same copynumber of fibers 
          
                /*logic_OpSurface_defaultAir[NofFiberscolumn][NofFibersrow] =
                    new G4LogicalBorderSurface("logic_OpSurface_defaultAir",
                                               CalorimeterPV, 
                                               physi_SiPM[column][row],
                                               OpSurfacedefault);*/
            }
        };
    };

    // Tubes with Cherenkov fibers and SiPM next to them
    //
    //G4VPhysicalVolume* physi_C_fiber[NofFiberscolumn][NofFibersrow];
  
    for(int column=0; column<NofFiberscolumn; column++){
        
        std::stringstream C_fiber_column;
        C_fiber_column.str("");
        C_fiber_column << column;
        for(int row=0; row<NofFibersrow; row++){
            
            std::stringstream C_fiber_row;
            C_fiber_row.str("");
            C_fiber_row << row;
            std::string C_name;
            std::string SiPM_name;
            C_name = "C_column_" + C_fiber_column.str() + "_row_" + C_fiber_row.str(); 
            SiPM_name = "C_SiPM"; 

            G4double C_x, C_y;
            G4ThreeVector vec_C_fiber;
            G4ThreeVector vec_SiPM;

            if(row%2 != 0){
                C_x = moduleX/2 - tuberadius - tuberadius - (tuberadius*2+2*tolerance)*column;
                C_y = -moduleY/2 + tuberadius + (sq3*tuberadius+2*tolerance)*row+tuberadius*(2.*sq3m1-1.);

                vec_C_fiber.setX(C_x);
                vec_C_fiber.setY(C_y);
                vec_C_fiber.setZ(0.);

                vec_SiPM.setX(C_x);
                vec_SiPM.setY(C_y);
                vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);

                copynumber = ((NofFibersrow/2)*column+row/2);
                        
                auto logic_C_fiber = constructcherfiber(tolerance,
                                                        tuberadius,
                                                        fiberZ,
                                                        absorberMaterial,
                                                        coreradius,
                                                        coreZ,
                                                        CherMaterial,
                                                        claddingradiusmin,
                                                        claddingradiusmax,
                                                        claddingZ,
                                                        CladCherMaterial);
                /*physi_C_fiber[column][row] =*/ new G4PVPlacement(0,
                                                         vec_C_fiber,
                                                         logic_C_fiber,
                                                         C_name,
                                                         moduleLV,
                                                         false,
                                                         copynumber);

                /*physi_SiPM[column][row] =*/ new G4PVPlacement(0,
                                                        vec_SiPM,
                                                        SiPMLV,
                                                        SiPM_name,
                                                        moduleequippedLV,
                                                        false,
                                                        copynumber); //same copynumber of fiber 

                /*logic_OpSurface_defaultAir[NofFiberscolumn][NofFibersrow] =
                    new G4LogicalBorderSurface("logic_OpSurface_defaultAir",
                                               CalorimeterPV, 
                                               physi_SiPM[column][row],
                                               OpSurfacedefault);*/
            }      
        };
    };




  G4VSolid* leakBoxS = new G4Box("leakbox", leakBoxX/2, leakBoxY/2, leakBoxZ/2);
  G4String sideString[4] = {"up", "right", "down", "left"};
  G4RotationMatrix siderotm = rotm;
  siderotm.rotateZ(90*deg);
  if(LeakageCounterIn)
  {
    for(int leakCounter = 0; leakCounter<NofLeakCounterLayers; leakCounter++)
    {
      //G4ThreeVector leakupPosition; leakupPosition.setX(fXshift); leakupPosition.setY(fYshift+caloY+leakBoxY+1.*cm); leakupPosition.setZ(-caloZ/2 - leakBoxZ/2 + leakBoxZ*leakCounter*2);
      G4ThreeVector leakupPosition; leakupPosition.setX(fXshift); leakupPosition.setY(fYshift+caloY+leakBoxY+distCaloBox); leakupPosition.setZ(-caloZ/2 - leakBoxZ/2 + caloZ/NofLeakCounterLayers*leakCounter*2);
      G4Transform3D transform_leakupBox = G4Transform3D(rotm, leakupPosition);
      G4LogicalVolume* leakupBoxLV = new G4LogicalVolume(leakBoxS, LeakCounterScinMaterial, "leakbox");
      G4VisAttributes* leakupBoxVisAtt = new G4VisAttributes(G4Colour(0.9,0,0.0));
      leakupBoxVisAtt->SetVisibility(true);
      leakupBoxVisAtt->SetForceWireframe(true);
      leakupBoxVisAtt->SetForceSolid(true);
      leakupBoxLV->SetVisAttributes(leakupBoxVisAtt);
      //new G4PVPlacement(transform_leakupBox, leakupBoxLV, "leakupBox", CalorimeterLV, false, 0, fCheckOverlaps);
      new G4PVPlacement(0, leakupPosition, leakupBoxLV, "leakbox", CalorimeterBoxLV, false, 4*leakCounter, fCheckOverlaps);

      G4ThreeVector leakrightPosition; leakrightPosition.setX(fXshift+caloX+leakBoxY+distCaloBox); leakrightPosition.setY(fYshift); leakrightPosition.setZ(-caloZ/2 - leakBoxZ/2 + caloZ/NofLeakCounterLayers*leakCounter*2);
      //G4Transform3D transform_leakrightBox = G4Transform3D(siderotm, leakrightPosition);
      G4Transform3D transform_leakrightBox = G4Transform3D(G4RotationMatrix(0., 0., 90*deg), leakrightPosition);
      G4LogicalVolume* leakrightBoxLV = new G4LogicalVolume(leakBoxS, LeakCounterScinMaterial, "leakbox");
      G4VisAttributes* leakrightBoxVisAtt = new G4VisAttributes(G4Colour(0.9,0,0.0));
      leakrightBoxVisAtt->SetVisibility(true);
      leakrightBoxVisAtt->SetForceWireframe(true);
      leakrightBoxVisAtt->SetForceSolid(true);
      leakrightBoxLV->SetVisAttributes(leakrightBoxVisAtt);
      new G4PVPlacement(transform_leakrightBox, leakrightBoxLV, "leakbox", CalorimeterBoxLV, false, 4*leakCounter+1, fCheckOverlaps);

      G4ThreeVector leakdownPosition; leakdownPosition.setX(fXshift); leakdownPosition.setY(fYshift-caloY-leakBoxY-distCaloBox); leakdownPosition.setZ(-caloZ/2 - leakBoxZ/2 + caloZ/NofLeakCounterLayers*leakCounter*2);
      G4Transform3D transform_leakdownBox = G4Transform3D(rotm, leakdownPosition);
      G4LogicalVolume* leakdownBoxLV = new G4LogicalVolume(leakBoxS, LeakCounterScinMaterial, "leakbox");
      G4VisAttributes* leakdownBoxVisAtt = new G4VisAttributes(G4Colour(0.9,0,0.0));
      leakdownBoxVisAtt->SetVisibility(true);
      leakdownBoxVisAtt->SetForceWireframe(true);
      leakdownBoxVisAtt->SetForceSolid(true);
      leakdownBoxLV->SetVisAttributes(leakdownBoxVisAtt);
      new G4PVPlacement(0, leakdownPosition, leakdownBoxLV, "leakbox", CalorimeterBoxLV, false, 4*leakCounter+2, fCheckOverlaps);

      G4ThreeVector leakleftPosition; leakleftPosition.setX(fXshift-caloX-leakBoxY-distCaloBox); leakleftPosition.setY(fYshift); leakleftPosition.setZ(-caloZ/2 - leakBoxZ/2 + caloZ/NofLeakCounterLayers*leakCounter*2);
      //G4Transform3D transform_leakleftBox = G4Transform3D(siderotm, leakleftPosition);
      G4Transform3D transform_leakleftBox = G4Transform3D(G4RotationMatrix(0., 0., 90*deg), leakleftPosition);
      G4LogicalVolume* leakleftBoxLV = new G4LogicalVolume(leakBoxS, LeakCounterScinMaterial, "leakbox");
      G4VisAttributes* leakleftBoxVisAtt = new G4VisAttributes(G4Colour(0.9,0,0.0));
      leakleftBoxVisAtt->SetVisibility(true);
      leakleftBoxVisAtt->SetForceWireframe(true);
      leakleftBoxVisAtt->SetForceSolid(true);
      leakleftBoxLV->SetVisAttributes(leakleftBoxVisAtt);
      new G4PVPlacement(transform_leakleftBox, leakleftBoxLV, "leakbox", CalorimeterBoxLV, false, 4*leakCounter+3, fCheckOverlaps);

    }

    G4ThreeVector tailCatcherPosition; tailCatcherPosition.setX(fXshift); tailCatcherPosition.setY(fYshift); tailCatcherPosition.setZ(caloZ + leakBoxY/2 + tailCatcherDist);
    G4RotationMatrix tailrotm;
    tailrotm.rotateX(90*deg);
    G4Transform3D transform_tailCatcherBox = G4Transform3D(tailrotm, tailCatcherPosition);
    G4LogicalVolume* tailCatcherBoxLV = new G4LogicalVolume(leakBoxS, LeakCounterScinMaterial, "leakbox");
    G4VisAttributes* tailCatcherBoxVisAtt = new G4VisAttributes(G4Colour(0.9,0,0.0));
    tailCatcherBoxVisAtt->SetVisibility(true);
    tailCatcherBoxVisAtt->SetForceWireframe(true);
    tailCatcherBoxVisAtt->SetForceSolid(true);
    tailCatcherBoxLV->SetVisAttributes(tailCatcherBoxVisAtt);
    //new G4PVPlacement(transform_tailCatcherBox, tailCatcherBoxLV, "tailCatcherBox", CalorimeterLV, false, 0, fCheckOverlaps);
    new G4PVPlacement(transform_tailCatcherBox, tailCatcherBoxLV, "leakbox", worldLV, false, 4*NofLeakCounterLayers, fCheckOverlaps);
  
  }


















    // Return physical world
    //
    return fWorldPV;

}

// Define constructscinfiber method()
//
G4LogicalVolume* DREMTubesDetectorConstruction::constructscinfiber(double tolerance,
    G4double tuberadius, G4double fiberZ, G4Material* absorberMaterial, 
    G4double coreradius, G4double coreZ, G4Material* ScinMaterial, 
    G4double claddingradiusmin, G4double claddingradiusmax, G4double claddingZ,
    G4Material* CherMaterial){
  

    auto invisibleAttr = new G4VisAttributes();
    invisibleAttr->SetVisibility(false);
   
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()m
    std::uniform_real_distribution<> dis(0.0, tolerance);
    double outradiussmear = dis(gen);
    tuberadius = tuberadius+outradiussmear;
  
    // Tube for scintillating fibers
    //
    G4Tubs* S_fiber = new G4Tubs("S_fiber", 0., tuberadius, fiberZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_S_fiber = new G4LogicalVolume(S_fiber,
                                                         absorberMaterial,
                                                         "S_fiber");
    //logic_S_fiber->SetVisAttributes(invisibleAttr);


    G4Tubs* Abs_S_fiber = new G4Tubs("Abs_Scin_fiber", claddingradiusmax, tuberadius, fiberZ/2,0.,2.*pi);

    G4LogicalVolume* logic_Abs_S_fiber = new G4LogicalVolume(Abs_S_fiber,
                                                             absorberMaterial,
                                                             "Abs_Scin_fiber");
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                               G4ThreeVector(0.,0.,0.),
                                               logic_Abs_S_fiber,
                                               "Abs_Scin_fiber",
                                               logic_S_fiber,
                                               false,
                                               0,
                                               fCheckOverlaps);

    G4Tubs* Core_S_fiber = new G4Tubs("Core_S_fiber", 0., 
                                      coreradius, coreZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Core_S_fiber = new G4LogicalVolume(Core_S_fiber,
                                                              ScinMaterial,
                                                              "Core_S_fiber");

    G4VisAttributes* ScincoreVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.8)); //blue
    ScincoreVisAtt->SetVisibility(true);
    ScincoreVisAtt->SetForceWireframe(true);
    ScincoreVisAtt->SetForceSolid(true);
    logic_Core_S_fiber->SetVisAttributes(ScincoreVisAtt);
    //logic_Core_S_fiber->SetVisAttributes(invisibleAttr);

    G4ThreeVector vec_Core_S;
    vec_Core_S.setX(0.);
    vec_Core_S.setY(0.);
    vec_Core_S.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     vec_Core_S,
                                                     logic_Core_S_fiber,
                                                     "Core_S_fiber",
                                                     logic_S_fiber,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);
 

    G4Tubs* Clad_S_fiber = new G4Tubs("Clad_S_fiber", claddingradiusmin, 
        claddingradiusmax, claddingZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Clad_S_fiber = new G4LogicalVolume(Clad_S_fiber,
                                                              CherMaterial,
                                                              "Clad_S_fiber");

    G4VisAttributes* ScincladVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    //light blue
    ScincladVisAtt->SetVisibility(true);
    ScincladVisAtt->SetForceWireframe(true);
    ScincladVisAtt->SetForceSolid(true);
    //logic_Clad_S_fiber->SetVisAttributes(ScincladVisAtt);
    logic_Clad_S_fiber->SetVisAttributes(invisibleAttr);


    G4ThreeVector vec_Clad_S;
    vec_Clad_S.setX(0.);
    vec_Clad_S.setY(0.);
    vec_Clad_S.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     vec_Clad_S,
                                                     logic_Clad_S_fiber,
                                                     "Clad_S_fiber",
                                                      logic_S_fiber,
                                                      false,
                                                      0,
                                                      fCheckOverlaps);


    G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(0.6,0.3,0.3)); //blue
    TubeVisAtt->SetVisibility(true);
    TubeVisAtt->SetForceWireframe(true);
    TubeVisAtt->SetForceSolid(true);
    logic_Abs_S_fiber->SetVisAttributes(TubeVisAtt);
    //logic_Abs_S_fiber->SetVisAttributes(invisibleAttr);
    
    return logic_S_fiber;

}

// Define constructcherfiber() method
//
G4LogicalVolume* DREMTubesDetectorConstruction::constructcherfiber(double tolerance,
    G4double tuberadius, G4double fiberZ, G4Material* absorberMaterial,
    G4double coreradius, G4double coreZ, G4Material* CherMaterial, 
    G4double claddingradiusmin, G4double claddingradiusmax, G4double claddingZ, 
    G4Material* CladCherMaterial){ 
 


    auto invisibleAttr = new G4VisAttributes();
    invisibleAttr->SetVisibility(false);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()m
    std::uniform_real_distribution<> dis(0.0, tolerance);
    double outradiussmear = dis(gen);
    tuberadius = tuberadius+outradiussmear;

    G4Tubs* C_fiber = new G4Tubs("C_fiber", 0., tuberadius, fiberZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_C_fiber = new G4LogicalVolume(C_fiber,
                                                         absorberMaterial,
                                                         "C_fiber");

    //logic_C_fiber->SetVisAttributes(invisibleAttr);

    G4Tubs* Abs_C_fiber = new G4Tubs("Abs_Cher_fiber", claddingradiusmax, tuberadius, fiberZ/2,0.,2.*pi);

    G4LogicalVolume* logic_Abs_C_fiber = new G4LogicalVolume(Abs_C_fiber,
                                                             absorberMaterial,
                                                             "Abs_Cher_fiber");
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     G4ThreeVector(0.,0.,0.),
                                                     logic_Abs_C_fiber,
                                                     "Abs_Cher_fiber",
                                                     logic_C_fiber,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);

    G4Tubs* Core_C_fiber = new G4Tubs("Core_C_fiber", 0., coreradius, coreZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Core_C_fiber = new G4LogicalVolume(Core_C_fiber,
                                                              CherMaterial,
                                                              "Core_C_fiber");

    G4VisAttributes* ChercoreVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
    ChercoreVisAtt->SetVisibility(true);
    ChercoreVisAtt->SetForceWireframe(true);
    ChercoreVisAtt->SetForceSolid(true);
    logic_Core_C_fiber->SetVisAttributes(ChercoreVisAtt);
    //logic_Core_C_fiber->SetVisAttributes(invisibleAttr);

    G4ThreeVector vec_Core_C;
    vec_Core_C.setX(0.);
    vec_Core_C.setY(0.);
    vec_Core_C.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                    vec_Core_C,
                                                    logic_Core_C_fiber,
                                                    "Core_C_fiber",
                                                    logic_C_fiber,
                                                    false,
                                                    0,
                                                    fCheckOverlaps);

    G4Tubs* Clad_C_fiber = new G4Tubs("Clad_C_fiber", claddingradiusmin,
        claddingradiusmax, claddingZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Clad_C_fiber = new G4LogicalVolume(Clad_C_fiber,
                                                              CladCherMaterial,
                                                              "Clad_C_fiber");

    G4VisAttributes* ChercladVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    //yellow 
    ChercladVisAtt->SetVisibility(true);
    ChercladVisAtt->SetForceWireframe(true);
    ChercladVisAtt->SetForceSolid(true);
    logic_Clad_C_fiber->SetVisAttributes(ChercladVisAtt);
    //logic_Clad_C_fiber->SetVisAttributes(G4VisAttributes::Invisible);
    logic_Clad_C_fiber->SetVisAttributes(invisibleAttr);

    G4ThreeVector vec_Clad_C;
    vec_Clad_C.setX(0.);
    vec_Clad_C.setY(0.);
    vec_Clad_C.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     vec_Clad_C,
                                                     logic_Clad_C_fiber,
                                                     "Clad_C_fiber",
                                                     logic_C_fiber,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);

    G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(0.6,0.3,0.3)); //blue
    TubeVisAtt->SetVisibility(true);
    TubeVisAtt->SetForceWireframe(true);
    TubeVisAtt->SetForceSolid(true);
    logic_Abs_C_fiber->SetVisAttributes(TubeVisAtt);
    //logic_Abs_C_fiber->SetVisAttributes(invisibleAttr);

    return logic_C_fiber;

}
//
//   method to define polygonal contour of module to extrude
//
std::vector<G4TwoVector> DREMTubesDetectorConstruction::calcmod(double radius, int nrow, int ncol) {
   G4double dxlr=radius;
   G4double dtubeY=sq3*radius;
   G4double moduleX = ncol*2.*radius+radius;
   G4double side=radius*2*sq3m1;
   G4double moduleY = (nrow-1.)*dtubeY+4.*radius*sq3m1;
   int nyp=2*nrow+1;
   int nxp=2*ncol-2;
   double yp[10000];
   double xp[10000];
   yp[0]=0.;
   xp[0]=0.;
   for(int i=1;i<nyp;i++){
     yp[i]=((i-1)%2+1)*0.5*side+yp[i-1];
     xp[i]=-((i+1)/2)%2*dxlr;
   }
   for(int i=nyp;i<(nyp+nxp);i++){
     int j=i-nyp;
     yp[i]=yp[nyp-1]+0.5*(i%2)*side;
     xp[i]=xp[nyp-1]+dxlr*(j+1);
   }
   for(int i=(nyp+nxp);i<(2*nyp+nxp);i++){
     int j=i-nyp-nxp;
     yp[i]=yp[nyp-j];
     xp[i]=xp[nyp+nxp-1]+dxlr+((j+1)/2)%2*dxlr;
   }
   for(int i=(2*nyp+nxp);i<(2*nyp+2*nxp);i++){
     int j=i-2*nyp-nxp;	
     yp[i]=0.5*side*(i%2);
     xp[i]=xp[nyp+nxp-j-1];
   }
//   xp[2*nyp+2*nxp]=0.;
//   yp[2*nyp+2*nxp]=0.;
   std::vector<G4TwoVector> polygon1(2*nyp+2*nxp);
   for(int i=0;i<(2*nyp+2*nxp);i++){
     polygon1[i].set(-xp[i]+moduleX/2.-radius,yp[i]-moduleY/2.);
   }
   return polygon1;
}
