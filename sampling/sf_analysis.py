#!/usr/bin/env python3
"""Sampling fraction vs beam energy from merged HidraSim ntuples.

For each input file (one per energy, hadd-merged) computes, per event,
    f_S   = EnergyScin / E0            (energy deposited in scintillating fibres)
    f_C   = EnergyCher / E0            (energy deposited in clear/Cherenkov fibres)
    f_S'  = EnergyScin / EnergyDepCalo (independent of longitudinal/lateral leakage)
    f_C'  = EnergyCher / EnergyDepCalo
    containment = EnergyDepCalo / E0
and prints mean, RMS and error on the mean. Writes a CSV and a PNG plot.

Usage: sf_analysis.py out/sf_E6.root out/sf_E40.root out/sf_E200.root [-o sf]
Needs uproot, numpy, matplotlib (pip install --user uproot matplotlib).
"""
import argparse
import csv

import numpy as np
import uproot

COLS = ["EnergyScin", "EnergyCher", "EnergyDepCalo", "PrimaryParticleEnergy"]


def stats(x):
    return x.mean(), x.std(ddof=1), x.std(ddof=1) / np.sqrt(len(x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--out", default="sf")
    a = ap.parse_args()

    rows = []
    for fn in a.files:
        t = uproot.open(fn)["DREMTubesout"]
        d = t.arrays(COLS, library="np")
        e0 = d["PrimaryParticleEnergy"]  # MeV
        ecalo = d["EnergyDepCalo"]
        q = {
            "fS": d["EnergyScin"] / e0,
            "fC": d["EnergyCher"] / e0,
            "fS_dep": d["EnergyScin"] / ecalo,
            "fC_dep": d["EnergyCher"] / ecalo,
            "containment": ecalo / e0,
        }
        row = {"file": fn, "E0_GeV": float(np.median(e0)) / 1000.0, "nevents": len(e0)}
        for k, v in q.items():
            m, s, em = stats(v)
            row[k], row[k + "_rms"], row[k + "_err"] = m, s, em
        rows.append(row)
    rows.sort(key=lambda r: r["E0_GeV"])

    keys = list(rows[0].keys())
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print(f"{'E0[GeV]':>8} {'N':>6} {'fS[%]':>16} {'fC[%]':>16} "
          f"{'fS/dep[%]':>16} {'fC/dep[%]':>16} {'contain[%]':>16}")
    for r in rows:
        line = f"{r['E0_GeV']:8.1f} {r['nevents']:6d}"
        for k in ["fS", "fC", "fS_dep", "fC_dep", "containment"]:
            line += f" {100*r[k]:8.4f}±{100*r[k+'_err']:6.4f}"
        print(line)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    e = np.array([r["E0_GeV"] for r in rows])
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, (k, lab) in zip(axs, [("fS", "scint"), ("fC", "Cherenkov"), ("containment", "E_dep/E0")]):
        for kk, mk, sfx in [(k, "o", "/E0"), (k + "_dep", "s", "/E_dep")] if k != "containment" else [(k, "o", "")]:
            y = np.array([r[kk] for r in rows]) * 100
            ey = np.array([r[kk + "_err"] for r in rows]) * 100
            ax.errorbar(e, y, yerr=ey, marker=mk, ls="-", label=lab + sfx)
        ax.set_xscale("log")
        ax.set_xlabel("E0 [GeV]")
        ax.set_ylabel("[%]")
        ax.legend()
        ax.grid(alpha=0.3)
    fig.suptitle("HidraSim: sampling fraction vs energy (e+, point-like beam)")
    fig.savefig(a.out + ".png", dpi=120)


if __name__ == "__main__":
    main()
