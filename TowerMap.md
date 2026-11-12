# Per-tower geometry map (CSV)

The simulation can dump a CSV file describing every **active tower** (module) of the
calorimeter, one row per tower, including whether the tower is read out by **SiPM** or
**PMT** and how its signal maps onto the output vectors. It is the tower-level companion of
the per-fiber map documented in [`FiberMap.md`](FiberMap.md).

Unlike the fiber map (SiPM only), the tower map contains **all** active towers.

## How to produce it

Add the UI command before `/run/initialize` in your macro:

```
/tbgeo/towermap towers.csv
/run/initialize
...
/run/beamOn 1
```

The file is written once, when the geometry is built (master thread), so a single event is
enough. Both maps can be requested in the same run (`/tbgeo/fibermap` and
`/tbgeo/towermap`). On completion the program prints:

```
WriteTowerMap: wrote <N> tower rows to '<path>' (<Ncol> columns x <Nrow> rows).
```

## Reference frame and conventions

- Coordinates are the **module centre** in the **calorimeter-surface frame with no rotation
  applied** (neither the platform tilt nor the `irot` 90° module rotation). Units are mm.
- `global_row` / `global_col` are **geometric ranks** of the tower centres, origin
  **bottom-left** (`global_row = 0` = smallest `y`, `global_col = 0` = smallest `x`), the
  same convention as the fiber map.

## Columns

| Column | Type | Description |
| --- | --- | --- |
| `tower_id` | int | Tower/module ID (`modflag` value = module copy number). |
| `readout` | string | `SiPM` if the tower is in `SiPMMod[]`, otherwise `PMT`. |
| `sipm_slot` | int | 0-based position of the tower in `SiPMMod[]` (the `SiPMTower` index), or `-1` for PMT towers. |
| `x` | float (mm) | Tower centre x (calorimeter frame, no rotation). |
| `y` | float (mm) | Tower centre y (calorimeter frame, no rotation). |
| `global_row` | int | Geometric rank of `y` over all towers; `0` = bottom. |
| `global_col` | int | Geometric rank of `x` over all towers; `0` = left. |
| `mod_grid_col` | int | Column of the tower in the calorimeter module grid (raw loop index). |
| `mod_grid_row` | int | Row of the tower in the calorimeter module grid (raw loop index). |
| `output_scin_vector` | string | Output vector for the scintillating signal: `VectorSignals` (SiPM) or `VecSPMT` (PMT). |
| `output_cher_vector` | string | Output vector for the Cherenkov signal: `VectorSignalsCher` (SiPM) or `VecCPMT` (PMT). |
| `output_index_base` | int | Base index of this tower in its output vectors (see below). |
| `n_fibers` | int | Total number of tubes in the tower, `NofFiberscolumn*NofFibersrow` (scintillating + Cherenkov). |
| `module_dx` | float (mm) | Module footprint width (x). |
| `module_dy` | float (mm) | Module footprint height (y). |
| `module_z` | float (mm) | Longitudinal length of the module/fibers (`moduleZ`). |

### `output_index_base` semantics

- **SiPM tower**: the fibers fill `VectorSignals` / `VectorSignalsCher` **per fiber**,
  starting at `output_index_base = sipm_slot * NoFibersTower`. Fiber `k` of the tower
  (`copynumber = k`) is at index `output_index_base + k`. See the per-fiber map for the
  individual entries.
- **PMT tower**: the whole tower is summed into a **single** entry
  `VecSPMT[output_index_base]` / `VecCPMT[output_index_base]`, where
  `output_index_base = tower_id`.

## Example

For the currently active geometry the file has `70` tower rows (62 PMT + 8 SiPM). The 8 SiPM
towers `{17,22,27,32,37,42,47,52}` form a vertical stack (`mod_grid_col = 2`) with
`sipm_slot` `0…7` and `output_index_base` `0, 512, … , 3584`.

```
tower_id,readout,sipm_slot,x,y,global_row,global_col,mod_grid_col,mod_grid_row,output_scin_vector,output_cher_vector,output_index_base,n_fibers,module_dx,module_dy,module_z
17,SiPM,0,0.0000,-97.0480,5,2,2,5,VectorSignals,VectorSignalsCher,0,1024,129.0000,28.3057,2500.0000
10,PMT,-1,256.0000,-124.7760,4,4,0,4,VecSPMT,VecCPMT,10,1024,129.0000,28.3057,2500.0000
```

## Visualization

`analysis/plot_towers.py` draws the towers as rectangles (x, y in mm), optionally
colored and/or labeled by any column:

```
python analysis/plot_towers.py towers.csv --color readout --label tower_id
python analysis/plot_towers.py towers.csv --color output_index_base --cmap plasma -o towers.png
```

### Note on the active geometry

In the currently compiled `modflag` block, tower id `1` is assigned to **two** physical
modules (grid rows 1 and 4) and id `11` is unused. The tower map lists both `tower_id = 1`
rows (different `y` / `mod_grid_row`); note that their PMT signals therefore add into the
same `VecSPMT[1]` / `VecCPMT[1]` entry. This is a property of the geometry configuration,
not of the map.
