/***************************************************/
/* Select prototype geometry and ancillaries setup */
/***************************************************/
#include "G4SystemOfUnits.hh"
// 
// Ancillaries
//
// Set to false to exclude the corresponding ancillary from simulation
const G4int NofLeakCounterLayers = 4;
const G4bool PreShowerIn = false;
const G4bool LeakageCounterIn = true;
const G4bool TruthLeakageIn = true;




// TB2021 
/*
    const G4int NofmodulesX = 3; 
    const G4int NofmodulesY = 3;
    const G4int modflag[9]={3,2,1,5,0,4,8,7,6};
    const G4int NoModulesSiPM=1;
    const G4int SiPMMod[1]={0};
    const G4int NofFiberscolumn = 16;
    const G4int NofFibersrow = 20;
    const G4int NoModulesActive=9;
    const G4double moduleZ = (1000.)*mm;
    const G4bool irot=false;
*/


// HiDRa, original 
    /*
    const G4int NofmodulesX = 24; 
    const G4int NofmodulesY = 5;
    const G4int modflag[120]={-1,-1,-1,-1,-1,-1,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-1,-1,-1,-1,-1,-1,-1,
                              -1,-1,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,-1,-1, 
                              30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53, 
                              -1,-1,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,-1,-1, 
                              -1,-1,-1,-1,-1,-1,-1,74,75,76,77,78,79,80,81,82,83,-1,-1,-1,-1,-1,-1,-1}; 
    const G4int NoModulesSiPM=10;
    const G4int SiPMMod[10]={37,38,39,40,41,42,43,44,45,46};
    const G4int NofFiberscolumn = 64;
    const G4int NofFibersrow = 16;
    const G4int NoModulesActive=84;
    const G4double moduleZ = (2500.)*mm;
    const G4bool irot=true;
    const G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;
    */


// TB 24  
/*
    const G4int NofmodulesX = 3;
    const G4int NofmodulesY = 12;
    // This one is DRAGO (if irot is true)
    const G4int modflag[36]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
    const G4int NoModulesSiPM=3;
    const G4int SiPMMod[NoModulesSiPM]={16, 19, 22};
    const G4int NofFiberscolumn = 64;
    const G4int NofFibersrow = 16;
    const G4int NoModulesActive=36;
    const G4bool irot=false;
    const G4double moduleZ = (2500.)*mm;
    const G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;
*/




/*
// 70 modules -> TB25
const G4int NofmodulesX = 5;
const G4int NofmodulesY = 18;
// 80 miniM -> But different geometry-> move one minimodules from row to corner, horizontal placement
const G4int modflag[90]={-1, -1, 0, -1, -1,
                         -1,  1, 2, 3, -1, 
                         -1,  4, 5,  6, -1, 
                         -1,  7, 8,  9, -1,
                         10,  11, 12, 13, 14, 
                         15,  16, 17, 18, 19, 
                         20,  21, 22, 23, 24,
                         25,  26, 27, 28, 29,
                         30,  31, 32, 33, 34, 
                         35,  36, 37, 38, 39, 
                         40,  41, 42, 43, 44, 
                         45,  46, 47, 48, 49, 
                         50,  51, 52, 53, 54,
                         55,  56, 57, 58, 59, 
                         -1,  60, 61, 62, -1, 
                         -1,  63, 64, 65, -1,
                         -1,  66, 67, 68, -1,
                         -1, -1,  69, -1, -1};  
const G4int NoModulesSiPM=8;
const G4int SiPMMod[8]={17,22,27,32,37,42,47,52};
const G4int NofFiberscolumn = 64;
const G4int NofFibersrow = 16;
const G4int NoModulesActive=70;
const G4double moduleZ = (2500.)*mm;
const G4bool irot=false;
const G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;
*/

// 80 modules -> TB26 (HiDRa geometry, SiPMs placed in two towers)
const G4int NofmodulesX = 5;
const G4int NofmodulesY = 20;
const G4int modflag[100]={-1, -1, 0, -1, -1,
                            -1,  1, 2, 3, -1, 
                            -1,  4, 5,  6, -1, 
                            -1, 7, 8,  9, -1,
                            10, 11, 12, 13, 14, 
                            15,  16, 17, 18, 19, 
                            20, 21, 22, 23, 24,
                            25, 26, 27, 28, 29,
                            30, 31, 32, 33, 34, 
                            35, 36, 37, 38, 39, 
                            40, 41, 42, 43, 44, 
                            45, 46, 47, 48, 49, 
                            50, 51, 52, 53, 54,
                            55, 56, 57, 58, 59, 
                            60, 61, 62, 63, 64, 
                            65, 66, 67, 68, 69, 
                            -1, 70, 71, 72, -1, 
                            -1, 73, 74, 75, -1,
                            -1, 76, 77, 78, -1,
                            -1, -1, 79, -1, -1};  
const G4int NoModulesSiPM=10;
//const G4int SiPMMod[10]={27,32,37,42,47,52,53, 48, 43, 38};
const G4int SiPMMod[10]={27,32,37,42,47,52, 26, 31, 36, 41};
const G4int NofFiberscolumn = 64;
const G4int NofFibersrow = 16;
const G4int NoModulesActive=80;
const G4double moduleZ = (2500.)*mm;
const G4bool irot=false;
const G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;






/*
// 80 modules -> Hidra Final geometry
const G4int NofmodulesX = 5;
const G4int NofmodulesY = 20;
const G4int modflag[100]={-1, -1, 0, -1, -1,
                            -1,  1, 2, 3, -1, 
                            -1,  4, 5,  6, -1, 
                            -1, 7, 8,  9, -1,
                            10, 1, 12, 13, 14, 
                            15,  16, 17, 18, 19, 
                            20, 21, 22, 23, 24,
                            25, 26, 27, 28, 29,
                            30, 31, 32, 33, 34, 
                            35, 36, 37, 38, 39, 
                            40, 41, 42, 43, 44, 
                            45, 46, 47, 48, 49, 
                            50, 51, 52, 53, 54,
                            55, 56, 57, 58, 59, 
                            60, 61, 62, 63, 64, 
                            65, 66, 67, 68, 69, 
                            -1, 70, 71, 72, -1, 
                            -1, 73, 74, 75, -1,
                            -1, 76, 77, 78, -1,
                            -1, -1, 79, -1, -1};  
const G4int NoModulesSiPM=10;
const G4int SiPMMod[10]={17,22,27,32,37,42,47,52,57,62};
const G4int NofFiberscolumn = 64;
const G4int NofFibersrow = 16;
const G4int NoModulesActive=80;
const G4double moduleZ = (2500.)*mm;
const G4bool irot=false;
const G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;
*/


