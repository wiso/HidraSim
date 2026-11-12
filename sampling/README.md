# Sampling fraction vs energy

Point-like e+ beam along +z, calorimeter at `/tbgeo/xshift = yshift = 0`
(its geometric centre on the beam axis). By default the beam line passes
through the centre of the calorimeter *front face*, which with the default
2.5°/2.5° tilts is at (x, y) ≈ (+54.5, −54.5) mm; `--target centre` aims at
the 3D centre of the calorimeter box instead (x = y = 0).

- `make_macro.py`: writes one job macro (energy, events, seed, tilts).
- `sf_analysis.py`: reads hadd-merged ntuples (one per energy) and prints /
  plots EnergyScin/E0, EnergyCher/E0, the same over EnergyDepCalo, and the
  containment EnergyDepCalo/E0.
