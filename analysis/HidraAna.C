//**************************************************
// \file HidraAna.C
// \brief:  analysis skeleton for HidraSim ntuples
// \author: Giacomo Polesello (INFN Pavia) 
//          Edoardo Proserpio (Uni Insubria)
//          Andrea Pareti (UniPV and INFN Pavia, andrea.pareti@cern.ch)
// \start date: May 3, 2022
//**************************************************
//
////usage: root -l -b -q 'HidraAna.C(energy,"filename")'
///  where energy is the energy of the beam, and filename
//   the name of the data ntuple
//
//   It produces an histogram file 
//   Hidra+"energy"+.root
//
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2F.h>
#include <iostream>
#include <array>
#include <stdint.h>
#include <string>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
// include file with geometry of module
#include "HidraGeo.h"


// Scintillator on even rows (start from 0), Cerenkov on uneven rows
// row ID input ranges from 0 to 31 
// need to know fiber type to distinguish global row
const unsigned int grouping = 8;

bool IsSiPMModule(int moduleID)
{
  for(int i=0; i<NoModulesSiPM; i++){
    if(SiPMMod[i] == moduleID) return true;
  }
  return false;
}

double ModuleCenterX(int moduleID, int modcol[])
{
  return dtubeX*NofFiberscolumn*((NofmodulesX-1.0)/2.0 - modcol[moduleID]);
}

double ModuleCenterY(int moduleID, int modrow[])
{
  return dtubeY*NofFibersrow*(modrow[moduleID] - (NofmodulesY-1.0)/2.0);
}

bool IsInsideModuleCell(double x, double y, int moduleID, int modcol[], int modrow[])
{
  const double pitchX = dtubeX*NofFiberscolumn;
  const double pitchY = dtubeY*NofFibersrow;
  const double centerX = ModuleCenterX(moduleID, modcol);
  const double centerY = ModuleCenterY(moduleID, modrow);

  return (std::abs(x - centerX) < pitchX/2.0) &&
         (std::abs(y - centerY) < pitchY/2.0);
}

bool IsInsideAnySiPMModuleCell(double x, double y, int modcol[], int modrow[])
{
  for(int i=0; i<NoModulesSiPM; i++){
    if(IsInsideModuleCell(x, y, SiPMMod[i], modcol, modrow)) return true;
  }
  return false;
}

std::string MakeEventDisplayName(const std::string& prefix, unsigned int entry)
{
  std::ostringstream os;
  os << prefix << "_evt" << std::setw(6) << std::setfill('0') << entry;
  return os.str();
}

TH2F* CreateHybridEventDisplay(const std::string& name,
                               const std::string& title,
                               double detectorExtentX,
                               double detectorExtentY)
{
  return new TH2F(name.c_str(), title.c_str(),
                  NofmodulesX*NofFiberscolumn/grouping, -detectorExtentX, detectorExtentX,
                  NofmodulesY*NofFibersrow, -detectorExtentY, detectorExtentY);
}

void FillTowerPatch(TH2F* hist, int moduleID, int modcol[], int modrow[], double content)
{
  if(!hist) return;

  const double pitchX = dtubeX*NofFiberscolumn;
  const double pitchY = dtubeY*NofFibersrow;
  const double towerXMin = ModuleCenterX(moduleID, modcol) - pitchX/2.0;
  const double towerXMax = ModuleCenterX(moduleID, modcol) + pitchX/2.0;
  const double towerYMin = ModuleCenterY(moduleID, modrow) - pitchY/2.0;
  const double towerYMax = ModuleCenterY(moduleID, modrow) + pitchY/2.0;

  const int firstXBin = std::max(1, hist->GetXaxis()->FindBin(towerXMin + 1.e-9));
  const int lastXBin = std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(towerXMax - 1.e-9));
  const int firstYBin = std::max(1, hist->GetYaxis()->FindBin(towerYMin + 1.e-9));
  const int lastYBin = std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(towerYMax - 1.e-9));

  for(int xbin=firstXBin; xbin<=lastXBin; xbin++){
    for(int ybin=firstYBin; ybin<=lastYBin; ybin++){
      const double x = hist->GetXaxis()->GetBinCenter(xbin);
      const double y = hist->GetYaxis()->GetBinCenter(ybin);
      if(IsInsideAnySiPMModuleCell(x, y, modcol, modrow)) continue;
      hist->SetBinContent(xbin, ybin, content);
    }
  }
}

void GetSiPMcoordinate(int TowID, int rowID, int colID_original, double &SiPM_X, double &SiPM_Y, std::string fiber, unsigned int grouping, int modcol_sipm, int modrow_sipm)
{

    // Calculate tower offsets based on the actual column and row positions of the SiPM module
    // Module spacing in detector: X = dtubeX*NofFiberscolumn, Y = dtubeY*NofFibersrow
    // Placement follows: m_x = dtubeX*NofFiberscolumn*((NofmodulesX-1)/2 - column)
    //                   m_y = dtubeY*NofFibersrow*(row - (NofmodulesY-1)/2)
    double TowerOffsetX = dtubeX*NofFiberscolumn*((NofmodulesX-1.0)/2.0 - modcol_sipm);
    double TowerOffsetY = dtubeY*NofFibersrow*(modrow_sipm - (NofmodulesY-1.0)/2.0);
    
    // after grouping, there are nfibercolumns/grouping channels
    unsigned int channel = static_cast<unsigned int>(colID_original/grouping);
    // colID should be such that the corresponding coordinate is in the middle of the channel,
    // between 0 and grouping-1 -> For 8 fibre grouping should be 3.5 (cast it to a double)
    double colID = (static_cast<double>(grouping)-1)/2 + channel*grouping;
    
    if(fiber == "S"){                                                         
        //SiPM_X = +moduleX/2 - tuberadius - (tuberadius*2)*colID;
        //SiPM_Y = -moduleY/2 + tuberadius + (sq3*tuberadius)*rowID+tuberadius*(2.*sq3m1-1.);
        SiPM_X = +moduleX/2 - tuberadius - (tuberadius*2)*colID;
        SiPM_Y = -moduleY/2 + tuberadius + (sq3*tuberadius)*rowID+tuberadius*(2.*sq3m1-1.);        
        //std::cout << "S X: " << SiPM_X << " Y: " << SiPM_Y << std::endl;
    }    

    if(fiber == "C"){                                                         
        // Cherenkov fibers are shifted by one additional tube radius in X compared to scintillating fibers
        SiPM_X = +moduleX/2 - (tuberadius*2) - (tuberadius*2)*colID ;
        SiPM_Y = -moduleY/2 + tuberadius + (sq3*tuberadius)*rowID+tuberadius*(2.*sq3m1-1.);
        //std::cout << "C X: " << SiPM_X << " Y: " << SiPM_Y << std::endl;
    }

    SiPM_X = TowerOffsetX + SiPM_X;
    SiPM_Y = TowerOffsetY + SiPM_Y;  
}





void HidraAna(double energy, const string intup, unsigned int EventDisplayEvery = 100){
//Open ntuples
  string infile = "../build/"+intup;
  std::cout<<"Using file: "<<infile<<std::endl;
  char cinfile[infile.size() + 1];
  strcpy(cinfile, infile.c_str());
  auto simfile = new TFile(cinfile, "READ");
  auto *simtree = (TTree*)simfile->Get( "DREMTubesout" );
  std::cout << "\n Test 0 ongoing... \n" << std::endl;

  std::ostringstream os;
  os << energy;
  std::string enstr = os.str();
  string outfile="hidra"+enstr+".root";
  TFile f(outfile.c_str(), "RECREATE");
  TDirectory* eventDisplayDir = f.mkdir("EventDisplays");
  std::cout << "\n Test 1 ongoing... \n" << std::endl;

//
//  build vectors with row and column position of each of
//  the MiniModules
//
  int modcol[NofmodulesX*NofmodulesY];
  int modrow[NofmodulesX*NofmodulesY];
  for(int i=0;i<NofmodulesX*NofmodulesY;i++){
    int row=i/NofmodulesX;
    int col=i%NofmodulesX;
    int imod=modflag[i];
    if(imod>=0){
      modcol[imod]=col;
      modrow[imod]=row;
    }
  }
  std::cout << "\n Test 2 ongoing... \n" << std::endl;
// book histograms  
  double bmin=energy-0.4*sqrt(energy)*10.;
  double bmax=energy+0.4*sqrt(energy)*10.;
  auto sciene = new TH1F("sciene", "sciene",100,bmin,bmax);
  auto cerene = new TH1F("cerene", "cerene",100,bmin,bmax);
  auto totene = new TH1F("totene", "totene",100,bmin,bmax);
  auto totenec = new TH1F("totenec", "totenec",100,bmin,bmax);
  auto totdep = new TH1F("totdep", "totdep",100,0.,bmax);
  auto leakene = new TH1F("leakene", "leakene",100,0.,0.1);
  auto chidist = new TH1F("chidist", "chidist",100,0.,1.);
  auto mapcalo  = new TH2F("mapcalo", "mapcalo",NofmodulesX,0.,NofmodulesX,NofmodulesY,0.,NofmodulesY);
  auto CaloCoordinatesMap = new TH2F("CaloCoordinatesMap", "Calo Coordinates Map; X [mm]; Y [mm]", 
                                     NofmodulesX, -dtubeX*NofFiberscolumn*(NofmodulesX/2.0), dtubeX*NofFiberscolumn*(NofmodulesX/2.0),
                                     NofmodulesY, -dtubeY*NofFibersrow*(NofmodulesY/2.0), dtubeY*NofFibersrow*(NofmodulesY/2.0));
  auto SipmMapS = new TH2F("SipmMapS", "SipmS; Col; Row", NofmodulesX*NofFiberscolumn, 0, NofmodulesX*NofFiberscolumn, NofmodulesY*NofFibersrow/2, 0, NofmodulesY*NofFibersrow);
  auto SipmMapC = new TH2F("SipmMapC", "SipmC; Col; Row", NofmodulesX*NofFiberscolumn, 0, NofmodulesX*NofFiberscolumn, NofmodulesY*NofFibersrow/2, 0, NofmodulesY*NofFibersrow);

  // Front view histograms: X-Y coordinates of the fiber hits as seen from +z looking towards detector
  double detectorExtentX = dtubeX * NofFiberscolumn * NofmodulesX / 2.0;
  double detectorExtentY = dtubeY * NofFibersrow * NofmodulesY / 2.0;

  auto SciSiPMCoordinates = new TH2F("SciSiPMCoordinates", "Sci SiPM Coordinates; X [mm]; Y[mm]",  NofmodulesX*NofFiberscolumn/grouping, -detectorExtentX, detectorExtentX, NofmodulesY*NofFibersrow/2, -detectorExtentY, detectorExtentY);
  auto CerSiPMCoordinates = new TH2F("CerSiPMCoordinates", "Cer SiPM Coordinates; X [mm]; Y[mm]",  NofmodulesX*NofFiberscolumn/grouping, -detectorExtentX, detectorExtentX, NofmodulesY*NofFibersrow/2, -detectorExtentY, detectorExtentY);

  int nentries=simtree->GetEntries();
  std::cout<<"Entries "<<nentries<<std::endl;

//Allocate branch pointers
  int pdg; simtree->SetBranchAddress( "PrimaryPDGID", &pdg );
  double venergy; simtree->SetBranchAddress( "PrimaryParticleEnergy", &venergy );
  double lenergy; simtree->SetBranchAddress( "EscapedEnergyl", &lenergy );
  double denergy; simtree->SetBranchAddress( "EscapedEnergyd", &denergy );
  double edep; simtree->SetBranchAddress( "EnergyDepCalo", &edep );
  double Stot; simtree->SetBranchAddress( "NofPMTScinDet", &Stot );
  double Ctot; simtree->SetBranchAddress( "NofPMTCherDet", &Ctot );
  double PSdep; simtree->SetBranchAddress( "PSEnergy", &PSdep );
  double beamX; simtree->SetBranchAddress( "PrimaryX", &beamX );
  double beamY; simtree->SetBranchAddress( "PrimaryY", &beamY );
  vector<double>* TowerE = NULL; 
  simtree->SetBranchAddress( "VecTowerE", &TowerE );
  vector<double>* SPMT = NULL; 
  simtree->SetBranchAddress( "VecSPMT", &SPMT );
  vector<double>* CPMT = NULL; 
  simtree->SetBranchAddress( "VecCPMT", &CPMT );
  vector<double>* SSiPM = NULL; 
  simtree->SetBranchAddress( "VectorSignals", &SSiPM );
  vector<double>* CSiPM = NULL; 
  simtree->SetBranchAddress( "VectorSignalsCher", &CSiPM );
// 
  double chi=0.38;   
  //double sciPheGeV=217.501;
  //double cerPheGeV=54.1621;
  //double sciPheGeV=125.501;
  //double cerPheGeV=29.6;  

  //const double sciPheGeV = 119.001; // tb24
  //const double cerPheGeV = 29.4;  // tb24
  //const double sciPheGeV = 178.501; // 10m
  //const double cerPheGeV = 43; // 10m
  const double sciPheGeV = 116.755;
  const double cerPheGeV = 29.654;

  double elcont=1.005;
  double picont=1.028;

  // total number of photoelectrons (not calibrated)
  double sciphe_raw = 0;
  double cerphe_raw = 0;
// Loop on events 
  for( unsigned int i=0; i<simtree->GetEntries(); i++){
    double ecalo=energy-lenergy/1000;
    simtree->GetEntry(i);
    const bool writeEventDisplay = (EventDisplayEvery > 0) && (((i + 1) % EventDisplayEvery) == 0);
    TH2F* eventDisplaySci = nullptr;
    TH2F* eventDisplayCer = nullptr;

    if(writeEventDisplay){
      std::ostringstream sciTitle;
      sciTitle << "Scintillation event display, entry " << i
               << "; X [mm]; Y [mm]";
      eventDisplaySci = CreateHybridEventDisplay(
        MakeEventDisplayName("EventDisplaySci", i),
        sciTitle.str(),
        detectorExtentX,
        detectorExtentY
      );

      std::ostringstream cerTitle;
      cerTitle << "Cherenkov event display, entry " << i
               << "; X [mm]; Y [mm]";
      eventDisplayCer = CreateHybridEventDisplay(
        MakeEventDisplayName("EventDisplayCer", i),
        cerTitle.str(),
        detectorExtentX,
        detectorExtentY
      );
    }

    double totsci=0.;
    double totcer=0.;
    double tottow=0.;



// Sum energy over all MiniModules
    for(unsigned int j=0; j<SPMT->size(); j++){
      const double sciTowerContent = SPMT->at(j)/sciPheGeV;
      const double cerTowerContent = CPMT->at(j)/cerPheGeV;
      sciphe_raw += SPMT->at(j);
      cerphe_raw += CPMT->at(j);
      totsci+=sciTowerContent;
      totcer+=cerTowerContent;
      tottow+=TowerE->at(j);
      mapcalo->Fill(modcol[j],modrow[j],TowerE->at(j)/1000/nentries);
      double towerX = dtubeX * NofFiberscolumn * ((NofmodulesX-1.0)/2.0 - modcol[j]);
      double towerY = dtubeY * NofFibersrow * (modrow[j] - (NofmodulesY-1.0)/2.0);
      CaloCoordinatesMap->Fill(towerX, towerY, TowerE->at(j)/1000/nentries);

      if(writeEventDisplay && !IsSiPMModule(j)){
        FillTowerPatch(eventDisplaySci, j, modcol, modrow, sciTowerContent);
        FillTowerPatch(eventDisplayCer, j, modcol, modrow, cerTowerContent);
      }
    }

   for(unsigned int N=0; N<SSiPM->size(); N++){        // Loop over SiPMs - S Fibers
      double content = SSiPM->at(N)/sciPheGeV;
      totsci+=content;
      sciphe_raw += SSiPM->at(N);
      unsigned int towID = static_cast<unsigned int>( N/(NofFiberscolumn*NofFibersrow/2) );
      unsigned int SiPMID = N%(NofFiberscolumn*NofFibersrow/2);
      unsigned int colID = static_cast<unsigned int>(SiPMID/(NofFibersrow/2));
      unsigned int rowID = 2*static_cast<unsigned int>(SiPMID%(NofFibersrow/2)); 
      // Get actual module ID from SiPM tower index
      int actual_mod_id = SiPMMod[towID];
      int modcol_sipm = modcol[actual_mod_id];
      int modrow_sipm = modrow[actual_mod_id];
      // Get coordinate
      double SiPM_X, SiPM_Y;
      GetSiPMcoordinate(towID, rowID, colID, SiPM_X, SiPM_Y, "S", grouping, modcol_sipm, modrow_sipm);
      //std::cout << towID << "\tSiPM ID: " << SiPMID << "\tcolID: " << colID << "\trowID: " << rowID << "\tX: " << SiPM_X << "\tY: " << SiPM_Y << std::endl;
      //std::cout << "Tower: " << towID << "\tModule row: " << modrow_sipm << "\tSiPMID: " << SiPMID << "\tColumn: " << colID << "\tRow: " << rowID << "\tTotal row: " << towID*NofFibersrow+rowID << std::endl;
      //SipmMapS->Fill( modcol_sipm*NofFiberscolumn + colID, modrow_sipm*NofFibersrow+rowID, content); 
      SipmMapS->Fill( modcol_sipm*NofFiberscolumn + colID, modrow_sipm*NofFibersrow+rowID, content); 
      SciSiPMCoordinates->Fill(SiPM_X, SiPM_Y, content);
      if(writeEventDisplay){
        eventDisplaySci->Fill(SiPM_X, SiPM_Y, content);
      }
      //std::cout << "S channel: " << "x: " << SiPM_X << "\t y: " << SiPM_Y << "\t content: " << content << std::endl;
    }


   for(unsigned int N=0; N<CSiPM->size(); N++){        // Loop over SiPMs - C Fibers
      double content = CSiPM->at(N)/cerPheGeV;
      totcer+=content;
      cerphe_raw += CSiPM->at(N);
      unsigned int towID = static_cast<unsigned int>( N/(NofFiberscolumn*NofFibersrow/2) );
      unsigned int SiPMID = N%(NofFiberscolumn*NofFibersrow/2);
      unsigned int colID = static_cast<unsigned int>(SiPMID/(NofFibersrow/2));
      unsigned int rowID = 2*static_cast<unsigned int>(SiPMID%(NofFibersrow/2)) + 1; // Cerenkov fibres on odd rows 
      // Get actual module ID from SiPM tower index
      int actual_mod_id = SiPMMod[towID];
      int modcol_sipm = modcol[actual_mod_id];
      int modrow_sipm = modrow[actual_mod_id];
      SipmMapC->Fill( modcol_sipm*NofFiberscolumn + colID, modrow_sipm*NofFibersrow+rowID, content); 
      double SiPM_X, SiPM_Y;
      GetSiPMcoordinate(towID, rowID, colID, SiPM_X, SiPM_Y, "C", grouping, modcol_sipm, modrow_sipm);
      CerSiPMCoordinates->Fill(SiPM_X, SiPM_Y, content);
      if(writeEventDisplay){
        eventDisplayCer->Fill(SiPM_X, SiPM_Y, content);
      }
      //std::cout << "C channel: " << "x: " << SiPM_X << "\t y: " << SiPM_Y << "\t content: " << content << std::endl;

    }

    if(writeEventDisplay){
      if(eventDisplayDir){
        eventDisplayDir->cd();
        eventDisplaySci->Write();
        eventDisplayCer->Write();
        f.cd();
      }
      delete eventDisplaySci;
      delete eventDisplayCer;
    }

    sciene->Fill(totsci);        
    cerene->Fill(totcer);        
    totene->Fill(elcont*0.5*(totsci+totcer));   
    totenec->Fill(picont*(totsci-chi*totcer)/(1-chi));   
    totdep->Fill(tottow/1000.);   
    leakene->Fill(lenergy/1000/energy);   
    chidist->Fill((totsci-ecalo)/(totcer-ecalo));    
    //std::cout << "totSci: " << totsci << "\t totCer: " << totcer << std::endl;
    //break;

  }
  std::cout << "Phe/GeV Sci: " << sciphe_raw/(venergy/1000)/nentries << "\t Phe/GeV Cer: " << cerphe_raw/(venergy/1000)/nentries << std::endl;
  
  /*
  totenec->Fit("gaus","Q","");
  TF1 *fit1 = totenec->GetFunction("gaus");
  double peak1=fit1->GetParameter(1);
  double epeak1=fit1->GetParError(1);
  double rms1=fit1->GetParameter(2);
  double erms1=fit1->GetParError(2);
  //cout << " # " << energy << " " << peak1 << " " << epeak1 << " " << rms1 << " " << erms1 << endl;
  */
  f.Write();
  //

}

//**************************************************
