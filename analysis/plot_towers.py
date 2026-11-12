#!/usr/bin/env python3
"""Plot the towers (modules) of the HiDRa calorimeter from a tower-map CSV.

The CSV is produced by the simulation with the macro command
``/tbgeo/towermap <path>`` (see TowerMap.md). Each tower is drawn as a rectangle
centred at its (x, y) position with size (module_dx, module_dy); axes in mm.

Optionally color and/or write a label inside each module by giving the name of a
CSV column.

Examples
--------
    python analysis/plot_towers.py towers.csv --label tower_id
    python analysis/plot_towers.py towers.csv --color readout --label tower_id
    python analysis/plot_towers.py towers.csv --color output_index_base --cmap plasma -o towers.png
"""
import argparse

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def color_patches(ax, pc, df, column, cmap):
    """Color a PatchCollection by ``column``; add colorbar/legend as needed."""
    if column is None:
        pc.set_facecolor("#e6e6e6")
        pc.set_edgecolor("0.2")
        pc.set_linewidth(0.6)
        ax.add_collection(pc)
        return

    series = df[column]
    if pd.api.types.is_numeric_dtype(series):
        pc.set_array(series.to_numpy())
        pc.set_cmap(cmap or "viridis")
        pc.set_edgecolor("0.2")
        pc.set_linewidth(0.6)
        ax.add_collection(pc)
        cbar = ax.figure.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(column)
    else:
        cats = list(pd.unique(series.astype(str)))
        qual = plt.get_cmap(cmap or "tab10")
        lut = {c: qual(i % qual.N) for i, c in enumerate(cats)}
        pc.set_facecolor([lut[str(v)] for v in series])
        pc.set_edgecolor("0.2")
        pc.set_linewidth(0.6)
        ax.add_collection(pc)
        handles = [
            plt.Line2D([0], [0], marker="s", linestyle="", markeredgecolor="0.2",
                       markerfacecolor=lut[c], label=c)
            for c in cats
        ]
        ax.legend(handles=handles, title=column, loc="upper right", fontsize=8)


COLUMNS = """\
Columns available for --color / --label (from /tbgeo/towermap CSV):
  tower_id            tower/module id (modflag)
  readout             SiPM or PMT
  sipm_slot           position in SiPMMod[] (0-based), or -1 for PMT
  x, y                tower centre in mm (calorimeter frame, no rotation)
  global_row          geometric row rank over towers (0 = bottom)
  global_col          geometric column rank over towers (0 = left)
  mod_grid_col        module column in the calorimeter grid
  mod_grid_row        module row in the calorimeter grid
  output_scin_vector  VectorSignals (SiPM) or VecSPMT (PMT)
  output_cher_vector  VectorSignalsCher (SiPM) or VecCPMT (PMT)
  output_index_base   base index of the tower in its output vectors
  n_fibers            number of tubes in the tower
  module_dx           module width in mm
  module_dy           module height in mm
  module_z            module length in mm
(Any other numeric or categorical column present in the CSV also works.)"""


def main():
    p = argparse.ArgumentParser(
        description="Plot towers as rectangles (x, y in mm).",
        epilog=COLUMNS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("csv", help="tower-map CSV (from /tbgeo/towermap)")
    p.add_argument("--color", help="column name to color the modules by")
    p.add_argument("--label", help="column name to write inside each module")
    p.add_argument("--cmap", help="matplotlib colormap name")
    p.add_argument("--fontsize", type=float, default=6.0, help="label font size")
    p.add_argument("--figsize", type=float, nargs=2, default=(8.0, 10.0),
                   metavar=("W", "H"), help="figure size in inches")
    p.add_argument("-o", "--output", help="save figure to this file instead of showing")
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    for col in ("x", "y", "module_dx", "module_dy"):
        if col not in df.columns:
            p.error(f"column '{col}' not found in {args.csv}")
    for col in (args.color, args.label):
        if col is not None and col not in df.columns:
            p.error(f"column '{col}' not found in {args.csv}")

    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    patches = [
        Rectangle((x - dx / 2.0, y - dy / 2.0), dx, dy)
        for x, y, dx, dy in zip(df.x, df.y, df.module_dx, df.module_dy)
    ]
    color_patches(ax, PatchCollection(patches, match_original=False),
                  df, args.color, args.cmap)

    if args.label:
        for x, y, v in zip(df.x, df.y, df[args.label]):
            ax.text(x, y, str(v), ha="center", va="center", fontsize=args.fontsize)

    xpad = df.module_dx.max()
    ypad = df.module_dy.max()
    ax.set_xlim((df.x - df.module_dx / 2).min() - 0.2 * xpad,
                (df.x + df.module_dx / 2).max() + 0.2 * xpad)
    ax.set_ylim((df.y - df.module_dy / 2).min() - 0.2 * ypad,
                (df.y + df.module_dy / 2).max() + 0.2 * ypad)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    title = "Calorimeter towers"
    if args.color:
        title += f" — colored by {args.color}"
    ax.set_title(title)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"saved {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
