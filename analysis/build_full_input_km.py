"""
Build a full-model input folder from input_e1500_enz2 with new Vmax-Km slopes (enzyme and uptake, MMRT and Arrhenius).

Usage:
    python build_full_input_km.py [--src ../input_e1500_enz2] [--dst ../input_e1500_enz2_km25] [--enzyme 2.5] [--uptake 2.5] [--force]
Sets mmrt_Vmax_Km_slope and Vmax_Km to --enzyme, and Uptake_mmrt_Vmax_Km_slope and Uptake_Vmax_Km to --uptake (previous values 7 and 5).
Everything else is copied unchanged.
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
    ap.add_argument("--src", default=os.path.join(HERE, "..", "input_e1500_enz2"))
    ap.add_argument("--dst", default=os.path.join(HERE, "..", "input_e1500_enz2_km25"))
    ap.add_argument("--enzyme", type=float, default=2.5)
    ap.add_argument("--uptake", type=float, default=2.5)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    src, dst = os.path.abspath(a.src), os.path.abspath(a.dst)
    if os.path.exists(dst):
        if not a.force:
            sys.exit(f"{dst} exists; use --force to rebuild it")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    set_params(os.path.join(dst, "parameters.csv"), {
        "mmrt_Vmax_Km_slope": a.enzyme, "Vmax_Km": a.enzyme,
        "Uptake_mmrt_Vmax_Km_slope": a.uptake, "Uptake_Vmax_Km": a.uptake})
    print("built", dst)


if __name__ == "__main__":
    main()
