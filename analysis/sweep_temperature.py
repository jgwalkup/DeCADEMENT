"""
Run the minimal scenario across a flat-temperature gradient (in-process, temp copies only),
optionally crossed with one swept parameter -- the temperature-axis counterpart to
sweep_minimal.py. Reuses sweep_minimal.run_one(); only the climate.csv Temp column changes
per point (Psi/FM/Date untouched), mirroring run_multi_climate_simulations.py's --static-temp
but applied in-place to an already-flat template.

Usage:
    # Temperature alone -- reference curve at the parameter baseline
    python sweep_temperature.py ../input_minimal --temps 4 8 12 16 20 24 28 32 36 40 \\
        --csv ../analysis_output/maint_thermal/summary.csv \\
        --timeseries ../analysis_output/maint_thermal/timeseries.csv

    # Temperature x parameter (2D sensitivity) -- one run per (temp, value) pair
    python sweep_temperature.py ../input_minimal --temps 4 12 20 28 36 \\
        --param Maint_E_HSP_Ea_Slope --values -30000 -15000 0 \\
        --csv out_summary.csv --timeseries out_timeseries.csv

PARAM (if given) may be a parameters.csv row name, or the stem of a NAME_min / NAME_max pair.
Every run: single pulse, side x side grid, fixed seed, flat climate at the given temperature.
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep_minimal import run_one  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("--temps", type=float, nargs="+", required=True, help="flat temperatures (deg C) to run")
    ap.add_argument("--param", default=None, help="optional second swept parameter (name or NAME_min/_max stem)")
    ap.add_argument("--values", type=float, nargs="+", default=None, help="values for --param (required if --param given)")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                    help="fixed override applied to every run (repeatable)")
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--csv", default=None, help="save the end-of-run summary table")
    ap.add_argument("--timeseries", default=None, help="save daily time series of every run (long format)")
    a = ap.parse_args()
    if a.param and not a.values:
        ap.error("--param requires --values")

    src = os.path.abspath(a.input_folder)
    fixed = {k: float(v) for k, v in (s.split("=") for s in a.set)}
    param_values = a.values if a.param else [None]

    print(f"# temperature sweep | input={src} | grid={a.side}x{a.side} | days={a.days} | seed={a.seed} "
          f"| fixed={fixed} | param={a.param}")
    rows, ts = [], []
    for t in a.temps:
        for v in param_values:
            changes = dict(fixed)
            if a.param:
                changes[a.param] = v
            summary, series = run_one(src, changes, a.side, a.days, a.seed, temp=t)
            row = {"Temp": t}
            if a.param:
                row[a.param] = v
            rows.append({**row, **summary})
            df = pd.DataFrame(series)
            df["Temp"] = t
            if a.param:
                df[a.param] = v
            ts.append(df)
            tag = f"Temp={t:g}" + (f" {a.param}={v:g}" if a.param else "")
            print(f"  done {tag}", file=sys.stderr, flush=True)

    index_cols = ["Temp"] + ([a.param] if a.param else [])
    out = pd.DataFrame(rows).set_index(index_cols)
    print(out.to_string(float_format=lambda x: f"{x:.4g}"))
    if a.csv:
        out.to_csv(a.csv)
    if a.timeseries:
        full = pd.concat(ts, ignore_index=True)
        ordered = index_cols + ["day"] + [c for c in full.columns if c not in index_cols + ["day"]]
        full[ordered].to_csv(a.timeseries, index=False)


if __name__ == "__main__":
    main()
