# Per-fiber geometry map (CSV)

The simulation can dump a CSV file describing every **SiPM-readout fiber** of the
calorimeter, one row per fiber. It is meant as a lookup table to connect an entry of the
output vectors (`VectorSignals` / `VectorSignalsCher`) to the physical fiber that produced
it: its type, its position, and its indices.

Only fibers belonging to SiPM towers (the towers listed in `SiPMMod[]` in
`include/DREMTubesGeoPar.hh`) are written. PMT-readout towers are **not** included.

## How to produce it

Add the UI command before `/run/initialize` in your macro:

```
/tbgeo/fibermap fibers.csv
/run/initialize
...
/run/beamOn 1
```

The file is written once, when the geometry is built (master thread), so a run with a
single event is enough. If the command is not given, no file is written. The path is
relative to the working directory. On completion the program prints:

```
WriteFiberMap: wrote <N> SiPM fiber rows to '<path>' (<Ncol> columns x <Nrow> rows).
```

## Reference frame and conventions

- Coordinates are given in the **calorimeter-surface frame with no rotation applied**:
  neither the test-beam platform tilt (`/tbgeo/vertrot`, `/tbgeo/horizrot`) nor the
  internal 90° module rotation (`irot`) affect the values. `x` is horizontal, `y` is
  vertical, both in **mm**.
- `global_row` / `global_col` are **geometric ranks**: all distinct `y` (resp. `x`)
  positions of the SiPM fibers are sorted, and each fiber gets the index of its position.
  Origin is **bottom-left**: `global_row = 0` is the smallest `y` (bottom),
  `global_col = 0` is the smallest `x` (left). Positions closer than 0.3 mm are treated as
  the same row/column.
- Because the fibers are hex-packed, scintillating (even rows) and Cherenkov (odd rows)
  fibers are staggered by half a tube in `x`, so they occupy **interleaved** columns. As a
  consequence, a given `(global_row, global_col)` cell holds at most one fiber, and the
  exact bottom-left *corner* may be empty (the fiber with the smallest `y` is generally not
  the one with the smallest `x`).

## Columns

| Column | Type | Description |
| --- | --- | --- |
| `output_index` | int | Index into the output vector for this fiber, exactly as filled by the stepping action: `SiPMTower * NoFibersTower + copynumber`. Range `0 .. NoFibersTower*NoModulesSiPM - 1`. |
| `output_vector` | string | Which output vector this fiber fills: `VectorSignals` (scintillating) or `VectorSignalsCher` (Cherenkov). |
| `fiber_type` | char | `S` = scintillating, `C` = Cherenkov. Equivalent to `output_vector`. |
| `x` | float (mm) | Global x in the calorimeter-surface frame (no rotation). |
| `y` | float (mm) | Global y in the calorimeter-surface frame (no rotation). |
| `global_row` | int | Geometric rank of `y` over all SiPM fibers; `0` = bottom. |
| `global_col` | int | Geometric rank of `x` over all SiPM fibers; `0` = left. |
| `local_row` | int | Row index within the module, `0` = bottom. Equal to the fiber row index (`0 .. NofFibersrow-1`). |
| `local_col` | int | Column index within the module, `0` = left. Equal to `NofFiberscolumn-1 - column`. |
| `tower_id` | int | Tower/module ID (`modflag` value), i.e. the module copy number. |
| `mod_grid_col` | int | Column of the module in the calorimeter module grid. |
| `mod_grid_row` | int | Row of the module in the calorimeter module grid. |
| `module_x` | float (mm) | x of the module centre (calorimeter frame). |
| `module_y` | float (mm) | y of the module centre (calorimeter frame). |
| `copynumber` | int | Fiber copy number within its module (= `SiPMID`): `(NofFibersrow/2)*column + row/2`, range `0 .. NoFibersTower-1`. |
| `pv_name` | string | Geant4 physical-volume name of the fiber, e.g. `S_column_12_row_4`. |

Notes:
- `output_index = SiPMTower * NoFibersTower + copynumber`, where `SiPMTower` is the
  position of `tower_id` in `SiPMMod[]` (0-based). The same `copynumber` value appears once
  per SiPM tower (with a different `output_index`).
- All geometry sizes (`NoFibersTower`, `SiPMMod`, `NofFiberscolumn`, `NofFibersrow`, …) are
  compile-time constants selected by the active block in `include/DREMTubesGeoPar.hh`. The
  CSV always reflects the block that was compiled in.

## Example

For the currently active geometry (8 SiPM towers `{17,22,27,32,37,42,47,52}` stacked
vertically, `NoFibersTower = 512`) the file has `8192` data rows and a `128 × 128` grid.
First data rows:

```
output_index,output_vector,fiber_type,x,y,global_row,global_col,local_row,local_col,tower_id,mod_grid_col,mod_grid_row,module_x,module_y,copynumber,pv_name
0,VectorSignals,S,63.5000,-110.0455,0,127,0,63,17,2,5,0.0000,-97.0480,0,S_column_0_row_0
1,VectorSignals,S,63.5000,-106.5795,2,127,2,63,17,2,5,0.0000,-97.0480,1,S_column_0_row_2
2,VectorSignals,S,63.5000,-103.1135,4,127,4,63,17,2,5,0.0000,-97.0480,2,S_column_0_row_4
```

## Visualization

`analysis/plot_fibers.py` draws the fibers as circles (x, y in mm), optionally colored
and/or labeled by any column:

```
python analysis/plot_fibers.py fibers.csv --color fiber_type
python analysis/plot_fibers.py fibers.csv --color output_index --cmap viridis -o fibers.png
```
