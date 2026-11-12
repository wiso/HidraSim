#!/usr/bin/env python3
"""Plot the SiPM fibers of the HiDRa calorimeter from a fiber-map CSV.

The CSV is produced by the simulation with the macro command
``/tbgeo/fibermap <path>`` (see FiberMap.md). Each fiber is drawn as a circle
at its (x, y) position, with the axes in millimetres.

Optionally color and/or write a label inside each fiber by giving the name of a
CSV column.

Examples
--------
    python analysis/plot_fibers.py fibers.csv
    python analysis/plot_fibers.py fibers.csv --color fiber_type
    python analysis/plot_fibers.py fibers.csv --color output_index --cmap viridis
    python analysis/plot_fibers.py fibers.csv --label copynumber -o fibers.png
"""
import argparse

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle


def color_patches(ax, pc, df, column, cmap):
    """Color a PatchCollection by ``column``; add colorbar/legend as needed."""
    if column is None:
        pc.set_facecolor("#7fb3ff")
        pc.set_edgecolor("0.3")
        pc.set_linewidth(0.2)
        ax.add_collection(pc)
        return

    series = df[column]
    if pd.api.types.is_numeric_dtype(series):
        pc.set_array(series.to_numpy())
        pc.set_cmap(cmap or "viridis")
        pc.set_edgecolor("none")
        ax.add_collection(pc)
        cbar = ax.figure.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(column)
    else:
        cats = list(pd.unique(series.astype(str)))
        qual = plt.get_cmap(cmap or "tab10")
        lut = {c: qual(i % qual.N) for i, c in enumerate(cats)}
        pc.set_facecolor([lut[str(v)] for v in series])
        pc.set_edgecolor("0.3")
        pc.set_linewidth(0.2)
        ax.add_collection(pc)
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", markeredgecolor="0.3",
                       markerfacecolor=lut[c], label=c)
            for c in cats
        ]
        ax.legend(handles=handles, title=column, loc="upper right", fontsize=8)


COLUMNS = """\
Columns available for --color / --label (from /tbgeo/fibermap CSV):
  output_index    index into VectorSignals/VectorSignalsCher for this fiber
  output_vector   VectorSignals (S) or VectorSignalsCher (C)
  fiber_type      S = scintillating, C = Cherenkov
  x, y            fiber position in mm (calorimeter frame, no rotation)
  global_row      geometric row rank over SiPM fibers (0 = bottom)
  global_col      geometric column rank over SiPM fibers (0 = left)
  local_row       row within the module (0 = bottom)
  local_col       column within the module (0 = left)
  tower_id        tower/module id (modflag)
  mod_grid_col    module column in the calorimeter grid
  mod_grid_row    module row in the calorimeter grid
  module_x        module centre x in mm
  module_y        module centre y in mm
  copynumber      fiber copy number within the module (= SiPMID)
  pv_name         Geant4 physical-volume name, e.g. S_column_12_row_4
(Any other numeric or categorical column present in the CSV also works.)"""


def main():
    p = argparse.ArgumentParser(
        description="Plot SiPM fibers as circles (x, y in mm).",
        epilog=COLUMNS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("csv", help="fiber-map CSV (from /tbgeo/fibermap)")
    p.add_argument("--tower", type=int, nargs="+", metavar="ID",
                   help="only plot fibers of the given tower_id(s)")
    p.add_argument("--color", help="column name to color the fibers by")
    p.add_argument("--label", help="column name to write inside each fiber")
    p.add_argument("--radius", type=float, default=1.0,
                   help="fiber circle radius in mm (default: 1.0)")
    p.add_argument("--cmap", help="matplotlib colormap name")
    p.add_argument("--fontsize", type=float, default=4.0, help="label font size")
    p.add_argument("--figsize", type=float, nargs=2, default=(8.0, 8.0),
                   metavar=("W", "H"), help="figure size in inches")
    p.add_argument("-o", "--output", help="save figure to this file instead of showing")
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    for col in ("x", "y"):
        if col not in df.columns:
            p.error(f"column '{col}' not found in {args.csv}")
    for col in (args.color, args.label):
        if col is not None and col not in df.columns:
            p.error(f"column '{col}' not found in {args.csv}")

    if args.tower is not None:
        if "tower_id" not in df.columns:
            p.error(f"column 'tower_id' not found in {args.csv}")
        df = df[df.tower_id.isin(args.tower)]
        if df.empty:
            available = sorted(int(t) for t in pd.read_csv(args.csv).tower_id.unique())
            p.error(f"no fibers for tower_id {args.tower}; available: {available}")

    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    patches = [Circle((x, y), args.radius) for x, y in zip(df.x, df.y)]
    color_patches(ax, PatchCollection(patches, match_original=False),
                  df, args.color, args.cmap)

    if args.label:
        if len(df) > 3000:
            print(f"warning: labeling {len(df)} fibers, this may be slow/dense")
        for x, y, v in zip(df.x, df.y, df[args.label]):
            ax.text(x, y, str(v), ha="center", va="center", fontsize=args.fontsize)

    r = args.radius
    ax.set_xlim(df.x.min() - 2 * r, df.x.max() + 2 * r)
    ax.set_ylim(df.y.min() - 2 * r, df.y.max() + 2 * r)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    title = "SiPM fibers"
    if args.tower:
        title += f" — tower {', '.join(str(t) for t in args.tower)}"
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
