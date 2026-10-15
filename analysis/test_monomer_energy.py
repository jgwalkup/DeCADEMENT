"""
Effect of the substrate monomer's energy yield (E_resp, E_ferm in monomer_energy.csv, row Mon_IdealSub)
on CUE and on the HSP effect, across temperature. For each (E_resp, E_ferm) variant this copies the
input folder, rewrites that row, and runs sweep_hsp_effect.py on the copy (HSP run + no-HSP control).

Usage:
    python test_monomer_energy.py ../input_minimal ../analysis_output/monomer_energy \\
        --variants 730,230 1000,230 1500,230 2200,230 730,700 730,2000 \\
        --set HSP_Induci_Prod=0.3 --set HSP_Consti_Prod=0.0003 --set Maint_E_HSP_Ea_Slope=-60000 \\
        --set HSP_prod_E_cost=300 --set Maint_E_Rate_Const=2 --temps 4 8 12 16 20 24 28 32 36
Writes <out_dir>/E<resp>_F<ferm>.csv per variant.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--variants", nargs="+", required=True, help="E_resp,E_ferm pairs")
    ap.add_argument("--row", default="Mon_IdealSub")
    ap.add_argument("--temps", nargs="+", default=["4", "8", "12", "16", "20", "24", "28", "32", "36"])
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--seed", default="0")
    a = ap.parse_args()

    src = os.path.abspath(a.input_folder)
    os.makedirs(a.out_dir, exist_ok=True)
    for v in a.variants:
        er, ef = v.split(",")
        tmp = tempfile.mkdtemp()
        work = os.path.join(tmp, "in")
        shutil.copytree(src, work)
        path = os.path.join(work, "monomer_energy.csv")
        lines = open(path).read().split("\n")
        out = []
        for l in lines:
            c = l.split(",")
            out.append(f"{c[0]},{er},{ef},{c[3]}" if c[0] == a.row else l)
        open(path, "w").write("\n".join(out))
        csv = os.path.join(a.out_dir, f"E{er}_F{ef}.csv")
        cmd = [sys.executable, os.path.join(HERE, "sweep_hsp_effect.py"), work, "--temps", *a.temps,
               "--seed", a.seed, "--csv", csv]
        for s in a.set:
            cmd += ["--set", s]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.rmtree(tmp, ignore_errors=True)
        print("done", v, flush=True)


if __name__ == "__main__":
    main()
