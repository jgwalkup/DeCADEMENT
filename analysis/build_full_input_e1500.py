"""
Build the full-model input folder `input_e1500` from `input/` with the parameters settled in the minimal-model analyses.

Usage:
    python build_full_input_e1500.py [--src ../input] [--dst ../input_e1500] [--force]
        [--enz-mult M] [--enz-init MIN MAX] [--taxa-per-box X]
    Adopted full-model setup (input_e1500_enz2):
        python build_full_input_e1500.py --dst ../input_e1500_enz2 --enz-mult 2 --enz-init 0.004 0.005 --taxa-per-box 0.05

Changes (everything else is copied unchanged):
  1. Energy yields: every E_resp and E_ferm in monomer_energy.csv is multiplied by 2 and rounded to the nearest ten
     (a non-zero yield is never rounded to 0: E_ferm 2 -> 4). E_value is left unchanged.
  2. Maintenance / HSP block of the minimal model: Maint_E_Rate_Const 2.5, Maint_E_HSP_Ea_Intercept 45000,
     Maint_E_HSP_Ea_Slope -90000, HSP_prod_E_cost 300; HSP_Consti_Prod 0.0015 and HSP_Induci_Prod 0.05 pegged (min = max).
  3. Enzyme and transporter production: per-gene rate ranges rescaled (relative range width kept) so the community-mean
     per-cell totals (rate x number of genes) equal the minimal model's:
        inducible enzymes    0.002  : Enz_Prod      E[rate] x 20 genes (Enz_per_taxon uniform 0-40)  -> E[rate] = 1.0e-4
        constitutive enzymes 0.0002 : Constit_Prod  E[rate] x 20 genes                              -> E[rate] = 1.0e-5
        transporter cost     0.225  : Uptake_C_cost E[cost] x 12.29 genes (measured over 6 communities) -> E[cost] = 0.01831
"""
import argparse
import os
import re
import shutil
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sweep_minimal import set_params  # noqa: E402

PARAMS = {
    "Maint_E_Rate_Const": 2.5,
    "Maint_E_HSP_Ea_Intercept": 45000,
    "Maint_E_HSP_Ea_Slope": -90000,
    "HSP_prod_E_cost": 300,
    "HSP_Consti_Prod_min": 0.0015,
    "HSP_Consti_Prod_max": 0.0015,
    "HSP_Induci_Prod_min": 0.05,
    "HSP_Induci_Prod_max": 0.05,
    "Enz_Prod_min": 4.62e-05,
    "Enz_Prod_max": 1.54e-04,
    "Constit_Prod_min": 5e-06,
    "Constit_Prod_max": 1.5e-05,
    "Uptake_C_cost_min": 0.01221,
    "Uptake_C_cost_max": 0.02441,
}


def scale_yield(v):
    x = 2 * v
    r = int(round(x / 10.0)) * 10
    return r if r > 0 else int(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=os.path.join(HERE, "..", "input"))
    ap.add_argument("--dst", default=os.path.join(HERE, "..", "input_e1500"))
    ap.add_argument("--force", action="store_true", help="replace an existing destination folder")
    ap.add_argument("--enz-mult", type=float, default=1.0,
                    help="multiplier applied to the matched Enz_Prod and Constit_Prod ranges (calibration, default 1)")
    ap.add_argument("--enz-init", type=float, nargs=2, default=None, metavar=("MIN", "MAX"),
                    help="initial enzyme pool per enzyme (Enz_min, Enz_max); default: unchanged (0.04 0.05)")
    ap.add_argument("--taxa-per-box", type=float, default=None, help="taxa_per_box in runtime.txt; default: unchanged (0.20)")
    a = ap.parse_args()
    src, dst = os.path.abspath(a.src), os.path.abspath(a.dst)
    if os.path.exists(dst):
        if not a.force:
            sys.exit(f"{dst} exists; use --force to rebuild it")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    params = dict(PARAMS)
    for k in ("Enz_Prod_min", "Enz_Prod_max", "Constit_Prod_min", "Constit_Prod_max"):
        params[k] = float(f"{PARAMS[k] * a.enz_mult:.6g}")
    if a.enz_init:
        params["Enz_min"], params["Enz_max"] = a.enz_init
    set_params(os.path.join(dst, "parameters.csv"), params)
    if a.taxa_per_box is not None:
        rp = os.path.join(dst, "runtime.txt")
        txt = open(rp).read()
        new = re.sub(r"taxa_per_box\t[0-9.]+", f"taxa_per_box\t{a.taxa_per_box:g}", txt)
        open(rp, "w").write(new)
    path = os.path.join(dst, "monomer_energy.csv")
    m = pd.read_csv(path, index_col=0)
    for c in ("E_resp", "E_ferm"):
        m[c] = m[c].map(scale_yield)
    m.to_csv(path)
    print(m.to_string())


if __name__ == "__main__":
    main()
