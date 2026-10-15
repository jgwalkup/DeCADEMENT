"""
Fungal share of the reseeded community across pulses, from taxa.csv of full_climate_osmo_range.py output.

Usage: python fungal_share.py OUT_DIR [--side 50] [--input ../input_e1500_enz2_km25_osmorange2]
Re-initialises each seed's community (same random stream as the run) to get which taxa are fungal (fb = 1), then combines the
per-taxon cell counts on day 1 and day 365 with the initial cell mass (bacteria 1, fungi 25 C).
Prints, by pulse: fungal share of cell-taxa and of initial biomass on day 1, and of cell-taxa and mean cells on day 365.
"""
import argparse
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from grid import Grid  # noqa: E402
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_runtime  # noqa: E402


def fungal_flags(src, seed, side):
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=1, end_time=5))
    cwd = os.getcwd()
    os.chdir(w)
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        eco = Grid(rt, initialize_data(rt))
        T = int(rt.loc["n_taxa", 1])
        fb = np.asarray(eco.fb[:T])
        names = list(pd.unique(eco.Microbes.index))
        c0 = np.asarray(eco.Microbes_init["C"]).reshape(side * side, T)[0]
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame({"taxon": names, "fungal": fb, "c0": c0})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--side", type=int, default=50)
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange2"))
    a = ap.parse_args()
    tx = pd.read_csv(os.path.join(a.out_dir, "taxa.csv"))
    flags = pd.concat([fungal_flags(os.path.abspath(a.input), s, a.side).assign(seed=s) for s in sorted(tx["seed"].unique())])
    tx = tx.merge(flags, on=["seed", "taxon"])
    G = a.side ** 2
    tx["C_day1"] = tx["cells_day1"] * tx["c0"]
    g = tx.groupby(["pulse", "seed"])
    out = pd.DataFrame({
        "rows_per_cell_day1": g["cells_day1"].sum() / G,
        "init_biomass": g["C_day1"].sum() / G,
        "fungal_share_rows_day1": g.apply(lambda x: x.loc[x.fungal == 1, "cells_day1"].sum() / x["cells_day1"].sum(), include_groups=False),
        "fungal_share_C_day1": g.apply(lambda x: x.loc[x.fungal == 1, "C_day1"].sum() / x["C_day1"].sum(), include_groups=False),
        "fungal_share_rows_day365": g.apply(lambda x: x.loc[x.fungal == 1, "cells_day365"].sum() / x["cells_day365"].sum(), include_groups=False),
        "fungal_taxa_present_day1": g.apply(lambda x: int(((x.fungal == 1) & (x.cells_day1 > 0)).sum()), include_groups=False),
    })
    print("mean over communities, by pulse:")
    print(out.groupby("pulse").mean().round(3).to_string())
    print("\ninitial biomass per cell by community and pulse:")
    print(out["init_biomass"].unstack("seed").round(1).to_string())
    out.to_csv(os.path.join(a.out_dir, "fungal_share.csv"))


if __name__ == "__main__":
    main()
