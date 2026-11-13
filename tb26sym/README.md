# TB26 vertical-symmetry study of the PMT-module profiles

Request: 20 GeV e+, TB26 geometry, platform horizontal rotation -2.5 deg, uniform
45 x 45 mm beam (+-22.5 mm in x and y, trigger acceptance), ~10k events per configuration:

| tag | beam centre on the front face | vertical tilt |
|---|---|---|
| `tb26sym_510c_vt0` | centre of PMT module 510 | 0 |
| `tb26sym_510c_vt1` | centre of PMT module 510 | 1 deg |
| `tb26sym_510b511_vt0` | boundary 510/511 | 0 |
| `tb26sym_510b511_vt1` | boundary 510/511 | 1 deg |

Physics list FTFP_BERT, build `~/HidraSim-tb26/build` on proof (Geant4 11.2.2, LCG_106b),
50 jobs x 200 events per configuration, registered in `/storage_tmp/atlas/turra/hidra/sim/RUNS.md`.

## Conventions

- Geometry: TB26 block of `include/DREMTubesGeoPar.hh` (5 columns x 20 rows, 80 modules).
  In the calorimeter frame module (sim column c, sim row r; 0-based) is centred at
  x = -128 mm (c - 2), y = 27.728 mm (r - 9.5); the copy number is `modflag[5 r + c]`.
- Module names follow the TB data convention `<column><row>` (1-based, see
  `TBDataPreparation/2025_SPS/PMT/mappingPMT.hpp`): row = r + 1 (row 1 at the bottom),
  column 5 = sim column 0 (x = +256 mm). So 509, 510, 511, 512 = copy numbers 30, 35, 40, 45.
  The 510/511 boundary is the horizontal mid-plane of the calorimeter (y = 0).
- Rotations: `/tbgeo/horizrot -2.5 deg`, `/tbgeo/vertrot 1 deg` (code convention: positive
  vertrot lifts the back of the calorimeter). The platform rotates about the calorimeter
  centre, so the beam line is computed (`sampling/make_macro.py --target-x/--target-y`) to go
  through the requested point of the front face: target (256, -13.864) mm for 510, (256, 0) mm
  for the boundary.
- Per event the ntuple has `ImpactX/ImpactY` (straight-line impact of the primary on the front
  face, calorimeter frame, mm), `VertexX/VertexY` (gun plane, world frame), and `VecSPMT`,
  `VecCPMT` (p.e.) and `VecTowerE` (MeV) indexed by copy number.

`slim.py` writes a flat tree with the impact point and S, C, E of modules 509-512.
