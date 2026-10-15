"""Quick MMRT vs Arrhenius comparison of the settled set (HSP and no-HSP control) across temperature, 365 days, seed 0."""
import os
from multiprocessing import Pool

from sweep_minimal import run_one

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.abspath(os.path.join(HERE, "..", "input_minimal_e1500"))
A = os.path.abspath(os.path.join(HERE, "..", "analysis_output", "kinetics_arrhenius", "input_minimal_e1500_arrhenius"))
OFF = {"HSP_Induci_Prod": 0.0, "HSP_Consti_Prod": 0.0}
TEMPS = (4, 10, 15, 20, 28, 36, 40)


def job(a):
    src, T, ch = a
    s, _ = run_one(src, ch, 20, 365, 0, temp=T)
    return s


if __name__ == "__main__":
    jobs = [(src, T, ch) for src in (M, A) for T in TEMPS for ch in ({}, OFF)]
    with Pool(12) as p:
        res = p.map(job, jobs, chunksize=1)
    i = 0
    print("kinetics  Temp  arm     deg%  biomC   CUE   E:C  alive")
    for lab in ("MMRT", "Arrh"):
        for T in TEMPS:
            for arm in ("HSP", "ctrl"):
                s = res[i]
                i += 1
                print(f"{lab}    {T:4d}  {arm:5s} {s['deg_pct']:6.2f} {s['biomC']:6.2f} {s['CUE']:5.2f} {s['EC_end']:6.0f} {s['alive']:4d}")
