#!/usr/bin/env python3
"""Slim per-event ntuple for the TB26 vertical-symmetry study (PyROOT, RDataFrame).

Input: merged HidraSim ntuple(s) (tree DREMTubesout) produced with the TB26 geometry
(DREMTubesGeoPar.hh, 5 columns x 20 rows).  Output tree "tb26sym" with, per event:
    ImpactX, ImpactY  primary impact point on the front face, calorimeter frame [mm]
    VertexX, VertexY  primary vertex in the gun plane, world frame [mm]
    S_<name>, C_<name> scintillation / Cherenkov signal [p.e.] of PMT module <name>
    E_<name>          energy deposited in the module [MeV]
    S_tot, C_tot, E_tot  sums over all modules
for the TB-data module names 509, 510, 511, 512 (column 5, rows 9-12).

Module name <col><row> (TB data convention, 1-based) -> simulation copy number:
    sim row    = row - 1     (row 1 at the bottom, y < 0)
    sim column = 5 - col     (sim column 0 is at x = +256 mm, "column 5"; use --mirror-x for
                              column 5 = sim column 4, x = -256 mm)
    copy number = modflag[sim row * 5 + sim column]
Module centres (calorimeter frame): x = -128 mm * (sim column - 2), y = 27.728 mm * (sim row - 9.5).

Usage: slim.py out.root in1.root [in2.root ...] [--mirror-x]
"""
import argparse

MODFLAG_TB26 = [-1, -1, 0, -1, -1,
                -1, 1, 2, 3, -1,
                -1, 4, 5, 6, -1,
                -1, 7, 8, 9, -1,
                10, 11, 12, 13, 14,
                15, 16, 17, 18, 19,
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
                -1, -1, 79, -1, -1]


def copy_number(name, mirror_x=False):
    col, row = int(name[0]), int(name[1:])
    scol = (col - 1) if mirror_x else (5 - col)
    cpn = MODFLAG_TB26[(row - 1) * 5 + scol]
    if cpn < 0:
        raise ValueError(f"module {name} does not exist")
    return cpn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--modules", default="509,510,511,512")
    ap.add_argument("--mirror-x", action="store_true")
    a = ap.parse_args()
    import ROOT
    df = ROOT.RDataFrame("DREMTubesout", a.inputs)
    cols = ["ImpactX", "ImpactY", "VertexX", "VertexY"]
    for name in a.modules.split(","):
        i = copy_number(name, a.mirror_x)
        for v, b in (("S", "VecSPMT"), ("C", "VecCPMT"), ("E", "VecTowerE")):
            df = df.Define(f"{v}_{name}", f"{b}[{i}]")
            cols.append(f"{v}_{name}")
        print(f"module {name} -> sim copy number {i}")
    df = (df.Define("S_tot", "Sum(VecSPMT)").Define("C_tot", "Sum(VecCPMT)")
            .Define("E_tot", "Sum(VecTowerE)"))
    cols += ["S_tot", "C_tot", "E_tot", "EnergyDepCalo" if "EnergyDepCalo" in df.GetColumnNames() else "EnergyTot"]
    df.Snapshot("tb26sym", a.out, cols)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
