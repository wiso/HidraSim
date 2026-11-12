# DREMTubes
**A Geant4 simulation of the 2020 Dual-Readout em-sized tubes prototype beam tests.**

<figure>
<img src="./images/DREMTubes_movie.gif" alt="Trulli" style="width:100%">
<figcaption align="center"><b>Fig. - 10 GeV positron passing through the preshower and the dual-readout calorimeter (2 events).</b></figcaption>
</figure>

<br/><br/>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-description">Project description</a></li>
    <li><a href="#authors-and-contacts">Authors and contacts</a></li>
     <li>
      <a href="#documentation-and-results">Documentation and results</a>
      <ul>
        <li><a href="#selected-presentations">Selected presentations</a></li>
      </ul>
    </li>
    <li><a href="#available-datasets-and-analyses">Available datasets and analyses</a></li>
    <li><a href="#output-variables">Output variables</a></li>
    <li>
      <a href="#how-to">How to</a>
      <ul>
        <li><a href="#build-compile-and-execute-on-maclinux">Build, compile and execute on Mac/Linux</a></li>
        <li><a href="#build-compile-and-execute-on-lxplus">Build, compile and execute on lxplus</a></li>
        <li><a href="#submit-a-job-with-htcondor-on-lxplus">Submit a job with HTCondor on lxplus</a></li>
      </ul>
    </li>
     </li><li><a href="#my-quick-geant4-installation">My quick Geant4 installation</a></li>
  </ol>                                           
</details>

<!--Project desription-->
## Project description
The project targets a standalone Geant4 simulation of the Dual-Readout tubes-based calorimeter prototype beam tests. We plan to perform Geant4 regression testing, physics lists comparison and validation against test-beam data. 
- Start date: 7 July 2021.

<!--Authors and contacts-->
## Authors and contacts
- (CERN EP-SFT) Lorenzo Pezzotti (lorenzo.pezzotti@cern.ch), Alberto Ribon (Supervisor)
- (University of Pavia and INFN Pavia) Andrea Pareti (andrea.pareti@cern.ch), Gabriella Gaudio

<!--Documentation and results-->
## Documentation and results

### Selected presentations
- Dual-Readout Calorimetry Meeting 19/11/2021, **Results from the CERN TB Geant4 simulation** [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://indico.cern.ch/event/1097245/contributions/4619316/attachments/2347762/4003816/lopezzot_DRSW_17_11_2021.pdf)
- Dual-Readout Calorimetry Meeting 13/10/2021, **Status of 2021 Test Beam(s) SW** [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://indico.cern.ch/event/1086651/contributions/4569695/attachments/2327255/3964777/lopezzot_DR_SW_13_10_2021.pdf)
- Dual-Readout Calorimetry Meeting 21/7/2021, **DREMTubes: A Geant4 simulation of the DR tubes prototype 2021 beam tests** [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://indico.cern.ch/event/1061304/contributions/4460441/attachments/2285253/3883980/DR_lopezzot_21_7_2021.pdf)


<!-- Output variables: -->
## Output variables

Each run writes a ROOT file named `DREMTubesout_Run<runID>.root`. The output ntuple is named `DREMTubesout` and is filled once per event. Values are written in Geant4 internal units, so energies are in MeV and positions or distances are in mm unless stated otherwise. Signal variables are simulated detected photoelectron counts, abbreviated as p.e.

### Detector-response parameters

The parameters controlling the simulated detector response are not stored in the output ntuple. They are defined in the source code and can be changed before recompiling:

| Parameter group | Location | Description |
| --- | --- | --- |
| Fast scintillation and Cherenkov signal smearing | `src/DREMTubesSignalHelper.cc` | The fast readout model converts deposited energy or optical photons to p.e. with Poisson smearing. The scintillation signal uses `SmearSSignalPMT()` and `SmearSSignalSiPM()`, while the Cherenkov signal uses `SmearCSignalPMT()` and `SmearCSignalSiPM()`. |
| Light attenuation lengths | `include/DREMTubesSignalHelper.hh` | `fSAttenuationLength` and `fCAttenuationLength` set the scintillation and Cherenkov attenuation lengths used by `AttenuateSSignal()` and `AttenuateCSignal()`. |
| Birks correction in the fast signal | `src/DREMTubesSignalHelper.cc` | `ApplyBirks()` applies the Birks correction before scintillation smearing in the fast signal model. |
| Material and optical properties | `src/DREMTubesDetectorConstruction.cc` | Refractive indices and optical surfaces are used by Geant4 optical photon processes, in particular for Cherenkov photon production/transport and boundary handling. The scintillation yield, emission spectrum, material Birks constant, and SiPM optical-surface efficiency/reflectivity are defined there, but they do not set the fast readout p.e. variables written to the ntuple; those are controlled by `DREMTubesSignalHelper`. Full optical propagation for scintillation is not supported in this version. |

### Scalar event variables

| Variable | Units | Description |
| --- | --- | --- |
| `EnergyScin` | MeV | Total energy deposited in scintillating fibers by ionizing charged particles with non-zero step length. It is accumulated step-by-step in scintillating fiber volumes before converting the signal to p.e. |
| `EnergyCher` | MeV | Total energy deposited in Cherenkov fiber volumes. It is accumulated step-by-step in Cherenkov fiber volumes. |
| `NofPMTCherDet` | p.e. | Total Cherenkov signal detected in PMT-readout towers. It is the event sum of all entries in `VecCPMT`. |
| `NofPMTScinDet` | p.e. | Total scintillation signal detected in PMT-readout towers. It is the event sum of all entries in `VecSPMT`. |
| `EnergyDepCalo` | MeV | Energy deposited in the calorimeter. It is the sum of `VecTowerE` (steel tubes, fiber cladding and core of all towers) plus the energy deposited in the air of the calorimeter structure (gaps between the tubes inside the modules and the `Calorimeter`/`CalorimeterBox` envelopes; about 0.5 MeV for a 20 GeV electron). The world volume, preshower volumes, leakage counters (`leakbox`) and truth-leakage absorber volumes are excluded. Invisible energy and energy escaping the detector (see `EscapedEnergyl`, `EscapedEnergyd`) are not included. |
| `PrimaryParticleEnergy` | MeV | Kinetic energy of the primary particle, saved from the primary track at its first step. |
| `PrimaryPDGID` | dimensionless | PDG particle ID of the primary particle. |
| `EscapedEnergyl` | MeV | Sum of kinetic energies of tracks entering the lateral truth-leakage absorber volume `leakageabsorberl`; the track is killed after being counted. |
| `EscapedEnergyd` | MeV | Sum of kinetic energies of tracks entering the longitudinal/downstream truth-leakage absorber volume `leakageabsorberd`; the track is killed after being counted. |
| `PSEnergy` | MeV | Total energy deposited in the preshower scintillator and lead volumes (if included in the simulation inside include/DREMTubesGeoPar.hh). |
| `PrimaryX` | mm | Primary-particle x position saved from the primary track at its first step. |
| `PrimaryY` | mm | Primary-particle y position saved from the primary track at its first step. |
| `NofSiPMScinDet` | p.e. | Total scintillation signal detected in SiPM-readout fibers. It is the event sum of all entries in `VectorSignals`. |
| `NofSiPMCherDet` | p.e. | Total Cherenkov signal detected in SiPM-readout fibers. It is the event sum of all entries in `VectorSignalsCher`. |

### Vector event variables

The vector lengths depend on the geometry selected in `include/DREMTubesGeoPar.hh`. The parameters used by the output vectors are:

| Geometry parameter | Description |
| --- | --- |
| `NofmodulesX` | Number of module slots in the x direction of the calorimeter module grid. |
| `NofmodulesY` | Number of module slots in the y direction of the calorimeter module grid. |
| `modflag` | Map from grid slot to active tower ID. The array has `NofmodulesX*NofmodulesY` entries; values below zero leave a slot empty, while non-negative values are used as module copy numbers and tower IDs. |
| `NoModulesActive` | Number of active towers/modules in the selected geometry. This sets the length of tower-level vectors such as `VecTowerE`, `VecSPMT`, and `VecCPMT`. |
| `NoModulesSiPM` | Number of active modules read out with SiPMs. This sets the module factor in the SiPM fiber-vector length. |
| `SiPMMod` | List of tower IDs that are read out with SiPMs. The position of a tower ID in this list is used as `SiPMTower` when indexing `VectorSignals` and `VectorSignalsCher`. Towers not listed here are treated as PMT-readout towers for `VecSPMT` and `VecCPMT`. |
| `NofFiberscolumn` | Number of fiber/tube columns in a module. |
| `NofFibersrow` | Number of fiber/tube rows in a module. Even rows are scintillating fibers and odd rows are Cherenkov fibers in the current placement logic. |
| `NoFibersTower` | Number of scintillating or Cherenkov fibers of one type in a module, defined as `NofFiberscolumn*NofFibersrow/2`. This sets the number of SiPM channels per SiPM-readout module for each signal type. |
| `NofLeakCounterLayers` | Number of leakage-counter layers on each lateral side. The leakage-counter vector length is `4*NofLeakCounterLayers+1`, including the downstream tail-catcher counter. |
| `PreShowerIn` | Enables or disables construction of the preshower volumes. If disabled, `PSEnergy` remains zero. |
| `LeakageCounterIn` | Enables or disables construction of the `leakbox` leakage counters. If disabled, `VecLeakCounter` is still created but remains zero. |
| `TruthLeakageIn` | Enables or disables the truth-leakage absorber volumes used by `EscapedEnergyl` and `EscapedEnergyd`. If disabled, those variables remain zero. |
| `moduleZ` | Longitudinal length of each calorimeter module and of the fibers inside it. |
| `irot` | Selects the orientation used when placing the module grid; it swaps the module x/y footprint used for calorimeter sizing and module placement. |

| Variable | Length | Units | Description |
| --- | ---: | --- | --- |
| `VectorSignals` | `NoModulesSiPM*NoFibersTower` | p.e. | Scintillation signal per SiPM-readout fiber. For each ionizing charged-particle step in a scintillating fiber, the deposited energy is Birks-corrected, converted to p.e. using a Poissonian smearing term, attenuated with scintillation attenuation length, and added at index `SiPMTower*NoFibersTower + SiPMID`. |
| `VectorSignalsCher` | `NoModulesSiPM*NoFibersTower` | p.e. | Cherenkov signal per SiPM-readout fiber. For optical photons in Cherenkov fibers undergoing total internal reflection, the signal is sampled with a Poissonian smearing term, attenuated with Cherenkov attenuation length, and added at index `SiPMTower*NoFibersTower + SiPMID`. |
| `VecTowerE` | `NoModulesActive` | MeV | Energy deposited in the tubes of each active tower, indexed by tower ID. For every scintillating and Cherenkov tube it includes the steel absorber tube (`Abs_Scin_fiber`, `Abs_Cher_fiber`), where most of the energy is deposited, and the fiber cladding and core. The air between the tubes is not included (it is counted only in `EnergyDepCalo`). |
| `VecSPMT` | `NoModulesActive` | p.e. | Scintillation signal per PMT-readout tower. It uses the same Birks correction, Poisson smearing, and attenuation as `VectorSignals`, but is filled only for towers that are not mapped to SiPM readout. |
| `VecCPMT` | `NoModulesActive` | p.e. | Cherenkov signal per PMT-readout tower. It uses the same Poisson smearing and attenuation as `VectorSignalsCher`, but is filled only for towers that are not mapped to SiPM readout. |
| `VecLeakCounter` | `4*NofLeakCounterLayers+1` | MeV | Energy deposited in leakage-counter `leakbox` volumes, indexed by leakage-counter copy number. For each layer, entries are ordered from the upper counter and then clockwise around the calorimeter: up, right, down, left. The first `4*NofLeakCounterLayers` entries correspond to the lateral counters, grouped by layer, and the final entry corresponds to the downstream tail catcher. |



<!--How to:-->
## How to

### Rotate and shift calorimeter
Note: the test-beam simulated platform can be shifted in x and y directions as in the actual configuration. The platform can also rotate around its center (the horizontal rotation). The housing containing the calorimeter can the lifted up from its back side creating a spin around its front face (the vertical rotation). By default, such parameters are set to zero. They are configurable via the UI macro card before the run is initialized as:
   ```
   /tbgeo/xshift <> [<Unit>]
   /tbgeo/yshift <> [<Unit>]
   /tbgeo/horizrot <> [<Unit>]
   /tbgeo/vertrot <> [<Unit>]
   ```

### Build, compile and execute on Mac/Linux
1. git clone the repo
   ```sh
   git clone https://github.com/DRCalo/HidraSim.git
   ```
2. source Geant4 env
   ```sh
   source /relative_path_to/geant4-v11.3.1-install/bin/geant4.sh
   ```
3. cmake build directory and make (using geant4-v11.3.1)
   ```sh
   mkdir build; cd build/
   cmake -DGeant4_DIR=/absolute_path_to/geant4-v11.3.1-install/lib64/Geant4-11.3.1/ relative_path_to/DREMTubes/
   make
   ```
4. execute (example with DREMTubes_run.mac macro card, 2 thread, FTFP_BERT physics list and no optical propagation)
   ```sh
   ./DREMTubes -m DREMTubes_run.mac -t 2 -pl FTFP_BERT
   ```
Parser options
   * -m macro.mac: pass a Geant4 macro card (example -m DREMTubes_run.mac available in source directory and automatically copied in build directory) 
   * -t integer: pass number of threads for multi-thread execution (example -t 3, default t=2)
   * -pl Physics_List: select Geant4 physics list (example -pl FTFP_BERT)
   * -opt FullOptic: boolean variable to switch on (true) the optical photon propagation in fibers (example -opt true, default false) -> NOTE: Not available any longer

### Build, compile and execute on lxplus and machines with CVMFS (ALMA9)
1. git clone the repo
   ```sh
   git clone git clone https://github.com/DRCalo/HidraSim.git
   ```
2. cmake, build directory and make (using geant4-v11.3.1, check for gcc and cmake dependencies for other versions)
   ```sh
   mkdir build; cd build/
   source /cvmfs/sft.cern.ch/lcg/views/LCG_106b/x86_64-el9-gcc11-opt/setup.sh
   cmake3 -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/11.3/x86_64-el9-gcc11-optdeb-MT/lib64/Geant4-11.3.0/ ../
   make (-jN)
   ```
3. execute (example with DREMTubes_run.mac macro card, 2 threads and FTFP_BERT physics list)
   ```sh
   ./DREMTubes -m DREMTubes_run.mac -t 2 -pl FTFP_BERT
   ```

### Build, compile and execute on local machine using vscode devcontainer

1. git clone the repo
   ```sh
   git clone https://github.com/lopezzot/DREMTubes.git
   cd DREMTubes
   ```

2. configure `.env` file
   ```sh
   cp .devcontainer/.env.example .devcontainer/.env
   ```
   in this new file edit the `GEANT4_DATASETS_HOST_PATH` variable with the path on your local machine where you want to store the Geant4 datasets (example: `GEANT4_DATASETS_HOST_PATH=$HOME/geant4-datasets`).

3. open the folder with vscode and open the devcontainer
   ```sh
   code .
   ```
   and click on "Reopen in container" when prompted, or open the command palette and search for "Dev Containers: Reopen in Container". This will build the docker image and start the container.

At this point you have an environment set up with Geant4 and all the needed dependencies. You should follow instructions in the terminal to build and execute the code.

### Submit a job with HTCondor on lxplus -> To be tested after Geant4-11
1. git clone the repo
   ```sh
   git clone https://github.com/lopezzot/DREMTubes.git
   ```
2. prepare execution files (example with Geant4.10.07_p01, DREMTubes_run.mac, 2 threads, FTFP_BERT physics list)
    ```sh
    mkdir DREMTubes-build; cd DREMTubes-build
    mkdir error log output
    cp ../../DREMTubes/scripts/DREMTubes_lxplus_10.7.p01.sh .
    source DREMTubes_lxplus_10.7.p01.sh
    ```
3. prepare for HTCondor submission (example with Geant4.10.07_p01, DREMTubes_run.mac, 2 threads, FTFP_BERT physics list)
    ```sh
    cp ../../DREMTubes/scripts/DREMTubes_HTCondor_10.7.p01.sh .
    export MYHOME=`pwd`
    echo cd $MYHOME >> DREMTubes_HTCondor_10.7.p01.sh
    echo $MYHOME/DREMTubes -m $MYHOME/DREMTubes_run.mac -t 2 >> DREMTubes_HTCondor_10.7.p01.sh
    cp ../../DREMTubes/scripts/DREMTubes_HTCondor.sub .
    sed -i '1 i executable = DREMTubes_HTCondor_10.7.p01.sh' DREMTubes_HTCondor.sub
    ```
4. submit a job
   ```sh
   condor_submit DREMTubes_HTCondor.sub 
   ```
5. monitor the job
   ```sh
   condor_q
   ```
   or (for persistency)
   ```sh
   condor_wait -status log/*.log
   ```

<!--My quick Geant4 installation-->
## My quick Geant4 installation
Here is my standard Geant4 installation (example with Geant4.10.7.p01) starting from the unpacked geant4.10.07.tar.gz file under the example path "path/to".

1. create build directory alongside source files
      ```sh
   cd /path/to
   mkdir geant4.10.07-build
   cd geant4.10.07-build
   ```
2. link libraries with CMAKE (example with my favourite libraries)
   ```sh
   cmake -DCMAKE_INSTALL_PREFIX=/Users/lorenzo/myG4/geant4.10.07_p01-install \
   -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_BUILD_MULTITHREADED=ON \
   -DGEANT4_USE_GDML=ON ../geant4.10.07.p01
   ```
3. make it
   ```sh
   make -jN
   make install
   ```
