"""
Build the full-model input folder from the final minimal parameter set (Phase 2.1 of SENSITIVITY_PLAN.md).

Usage:
    python analysis/build_full_input_port.py [--src input_e1500_enz2_km25_osmorange2_bac] [--dst input_full_port_bac]
        [--litter-scale 1.245] [--enzyme-factor 4.5] [--trait-width 1.25] [--force]

Starts from the latest bacterial full-model folder (fb 0) and applies:

Energy yields (monomer_energy.csv)
  - the ten litter monomers (cellulose ... OrgP2): E_resp x litter-scale (default 1.245), E_ferm and E_value unchanged.
    With the default, the litter-C-weighted maximum energy per unit C taken up, sum_s C_s f_s (0.9 E_resp + 0.1 E_ferm) /
    sum_s C_s with routing f_s = E_value ** ((2 / log10 E_ferm) ** Ferm_Cond), equals that of the minimal substrate
    (IdealSub E_resp 1500, E_ferm 230, E_value 1: 1373); the e1500 table gives 1109.
  - Mon_DeadMic, Mon_DeadEnz and FermProd: the minimal values (700 / 100 / 0.4, 350 / 5 / 0.4, 690 / 2 / 0.5), so endogenous
    respiration and the rebalancing credit (both at the Mon_DeadMic E_resp) have the minimal model's yield.

Parameters (parameters.csv), as in the final minimal baseline (SENSITIVITY_SUMMARY sections 29-40)
  - Maint_E_Rate_Const 4, Enz_loss_rate 0.05, Vmax_Km 6 (Uptake_Vmax_Km 0.25 unchanged), HSP_prod_E_cost 600,
    Osmo_prod_E_cost 150, Input_PO4 0.005 (beta 1, psi slopes 0.1, R_endo / R_stoich / Mineral_Uptake_Passive 1 are
    already set in the source folder).
  - Input_NH4 0.1 (minimal: 0.05): the full-model taxon grows faster and runs out of NH4 at 0.05 (section 41c-41d).
  - Enzyme production normalised per taxon: NormalizeProd 1, so each taxon's total production is the per-gene rate x
    Enz_per_taxon_max whatever its number of enzyme genes. Constit_Prod and Enz_Prod are set to enzyme-factor (default
    4.5) x the minimal totals (0.00045 and 0.008) / Enz_per_taxon_max (40): 5.0625e-5 and 9e-4 for every taxon
    (min = max). The factor makes up for enzyme specificity: each full-model enzyme acts on one substrate, so a taxon's
    decay capacity per unit production is about 7x lower than the minimal taxon's (section 41). The trait ranges for the
    community are set in step 2.3.
  - Community trait ranges (step 2.3; narrow to start): each taxon draws its rates uniformly (Latin hypercube) from
    centre / w .. centre x w with w = trait-width (default 1.25): Constit_Prod and Enz_Prod around the values above,
    Osmo_Consti_Prod around 0.005 and Osmo_Induci_Prod around 0.05 (minimal values), HSP_Consti_Prod around 0.0015 and
    HSP_Induci_Prod around 0.05; Enz_per_taxon 30-40 genes of the 50 (no cheaters). --trait-width 1 gives single
    values (min = max) and leaves the source's gene range.
  - Single-taxon check at these settings (section 41d; ambient, MMRT): 91 % of the litter degraded, day 50 % 130,
    day 95 % 313, CUE 0.37, endogenous respiration 34 % of respiration.
  - Arrhenius pre-exponential constants as in the minimal (matched to MMRT at 12.5-20 C on the main substrate):
    Vmax0 5e7, Uptake_Vmax0 4.2e7 (min = max; the source drew 1e7-1e8 per enzyme x substrate).
  - Not changed: osmolyte production ranges (community trait variation), HSP production, uptake gene costs, initial
    enzyme pools, MMRT HT0 / ST0 / Cp ranges, substrates and their inputs, runtime.txt.
"""
import argparse
import os
import shutil
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sweep_minimal import set_params  # noqa: E402

LITTER = ["Mon_Cellulose", "Mon_Hemicellulose", "Mon_Starch", "Mon_Chitin", "Mon_Lignin", "Mon_Protein1",
          "Mon_Protein2", "Mon_Protein3", "Mon_OrgP1", "Mon_OrgP2"]
NECRO = {"FermProd": (690, 2, 0.5), "Mon_DeadMic": (700, 100, 0.4), "Mon_DeadEnz": (350, 5, 0.4)}
MINIMAL_TOTAL = {"Constit_Prod": 0.00045, "Enz_Prod": 0.008}
PARAMS = {"Maint_E_Rate_Const": 4.0, "Enz_loss_rate": 0.05, "Vmax_Km": 6.0, "HSP_prod_E_cost": 600.0,
          "Osmo_prod_E_cost": 150.0, "Input_NH4": 0.1, "Input_PO4": 0.005, "NormalizeProd": 1,
          "Vmax0": 5e7, "Uptake_Vmax0": 4.2e7}
FERM_COND = 0.1


def weighted_yield(me, subs):
    """Litter-C-weighted maximum energy per unit C taken up."""
    import numpy as np
    s = t = 0.0
    for sub, c in subs["C"].items():
        m = "Mon_" + sub
        if m not in LITTER:
            continue
        r = me.loc[m]
        f = r.E_value ** ((2 / np.log10(r.E_ferm)) ** FERM_COND)
        s += c * f * ((1 - FERM_COND) * r.E_resp + FERM_COND * r.E_ferm)
        t += c
    return s / t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange2_bac"))
    ap.add_argument("--dst", default=os.path.join(HERE, "..", "input_full_port_bac"))
    ap.add_argument("--litter-scale", type=float, default=1.245)
    ap.add_argument("--enzyme-factor", type=float, default=4.5)
    ap.add_argument("--trait-width", type=float, default=1.25)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    src, dst = os.path.abspath(a.src), os.path.abspath(a.dst)
    if os.path.exists(dst):
        if not a.force:
            sys.exit(f"{dst} exists; use --force to rebuild it")
        shutil.rmtree(dst)
    shutil.copytree(src, dst, symlinks=False, ignore=shutil.ignore_patterns(".DS_Store", ".RData", ".Rhistory"))

    # energy yields
    p_me = os.path.join(dst, "monomer_energy.csv")
    me = pd.read_csv(p_me, index_col=0).astype(float)
    subs = pd.read_csv(os.path.join(dst, "initial_substrates.csv"), index_col=0)
    before = weighted_yield(me, subs)
    me.loc[LITTER, "E_resp"] = (me.loc[LITTER, "E_resp"] * a.litter_scale).round(0)
    for m, v in NECRO.items():
        me.loc[m, ["E_resp", "E_ferm", "E_value"]] = v
    me.to_csv(p_me, float_format="%g")
    after = weighted_yield(me, subs)

    # parameters
    p_par = os.path.join(dst, "parameters.csv")
    par = pd.read_csv(p_par, header=None, index_col=0)[1]
    n_max = float(par["Enz_per_taxon_max"])
    changes = dict(PARAMS)
    for k, tot in MINIMAL_TOTAL.items():
        changes[k] = float(f"{a.enzyme_factor * tot / n_max:.6g}")
    if a.trait_width != 1:
        w = a.trait_width
        centres = {k: changes.pop(k) for k in MINIMAL_TOTAL}
        centres.update({"Osmo_Consti_Prod": 0.005, "Osmo_Induci_Prod": 0.05, "HSP_Consti_Prod": 0.0015,
                        "HSP_Induci_Prod": 0.05})
        for k, c in centres.items():
            changes[k + "_min"] = float(f"{c / w:.6g}")
            changes[k + "_max"] = float(f"{c * w:.6g}")
        changes["Enz_per_taxon_min"], changes["Enz_per_taxon_max"] = 30, 40
    set_params(p_par, changes)

    with open(os.path.join(dst, "README.md"), "w") as fh:
        fh.write(f"# {os.path.basename(dst)}\n\nBuilt by `analysis/build_full_input_port.py` from "
                 f"`{os.path.basename(src)}/` (litter scale {a.litter_scale}, enzyme factor {a.enzyme_factor}, trait width "
                 f"{a.trait_width}). See the script docstring for every "
                 f"change.\n\nLitter-C-weighted maximum energy per unit C taken up: {before:.0f} (source) -> "
                 f"{after:.0f} (minimal IdealSub: 1373).\n")
    print(f"built {dst}; litter-weighted energy per uptake C {before:.0f} -> {after:.0f}; set "
          + ", ".join(f"{k} {v:g}" for k, v in changes.items()))


if __name__ == "__main__":
    main()
