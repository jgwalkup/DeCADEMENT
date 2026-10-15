"""
Copy a full-model input folder and give the taxa a range of osmolyte production (so relative drought tolerance varies).

Usage:
    python build_full_input_osmo_range.py --src ../input_e1500_enz2_km25 --dst ../input_e1500_enz2_km25_osmorange
        [--induci 0.01 0.09] [--consti 0.005 0.005] [--force]
Sets Osmo_Induci_Prod_min/max and Osmo_Consti_Prod_min/max; each taxon draws its own rates (Latin hypercube, uniform), so its total
allocation is constitutive + inducible. Everything else, including Osmo_tol_lower and Osmo_tol_ref, is copied unchanged.
"""
import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sweep_minimal import set_params  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=os.path.join(HERE, "..", "input_e1500_enz2_km25"))
    ap.add_argument("--dst", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange"))
    ap.add_argument("--induci", type=float, nargs=2, default=[0.01, 0.09])
    ap.add_argument("--consti", type=float, nargs=2, default=[0.005, 0.005])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    src, dst = os.path.abspath(a.src), os.path.abspath(a.dst)
    if os.path.exists(dst):
        if not a.force:
            sys.exit(f"{dst} exists; use --force to rebuild it")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    set_params(os.path.join(dst, "parameters.csv"), {
        "Osmo_Induci_Prod_min": a.induci[0], "Osmo_Induci_Prod_max": a.induci[1],
        "Osmo_Consti_Prod_min": a.consti[0], "Osmo_Consti_Prod_max": a.consti[1]})
    print("built", dst)


if __name__ == "__main__":
    main()
