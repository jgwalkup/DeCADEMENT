"""
-------------------------------------------------------------------------------
         DEMENTpy--Decomposition Model of Enzymatic Traits in Python
                      Modified to calculate rates using MMRT
                              Jeth Walkup, PhD
              Department of Ecology and Evolutionary Biology
                       University of California Irvine
                  Emails: jethwalkup@gmail.com or jwalkup@uci.edu
-------------------------------------------------------------------------------
"""

import os
import sys, traceback
import time
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

import pandas as pd
import numpy  as np

from initialization import initialize_data
from grid import Grid
from output import Output
from utility import export


def main():
    
    
    print("""
    ---------------------------------------------------------------------------
         DEMENTpy (DEcomposition Model of Enzymatic Traits in Python)
                               MMRT Version 
               Department of Ecology and Evolutionary Biology
                     University of California Irvine
    ---------------------------------------------------------------------------       
    """)
    
     
    #...Obtain the command line arguments
    input_folder  = sys.argv[1]                                    # input folder name
    output_folder = sys.argv[2]                                    # output folder name
    outname       = sys.argv[3]                                    # output file name / run label
    climate_file  = sys.argv[4] if len(sys.argv) > 4 else 'climate.csv'  # optional climate file name

    #...Set up the working directory
    os.chdir('../'+input_folder)

    #...grow a seed of random number generator from the digits in outname
    digits = ''.join(c for c in outname if c.isdigit())
    np.random.seed(int(digits) % (2**32) if digits else 0)

    #...a few system constants
    runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
    output_format = str(runtime.loc['output_format', 1]).strip().lower()
    runtime = runtime.drop('output_format')
    runtime = runtime.astype(float)
    pulse      = int(runtime.loc['pulse',1])         # number of pulses
    cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
    # Export_format = str(runtime.loc['output_format', 1]).strip().lower() # output file format (CSV, NetCDF, Zarr)
    mic_reinit = True    # indicate reinitialization of microbial community
    
    #...Initialize data by calling the Function: Initialize_Data()
    data_initialization = initialize_data(runtime, climate_file)

    #...Prepare for output by creating an instance of the Output class
    Output_init = Output(runtime,data_initialization)

    #...Create an instance of the Grid class
    Ecosystem = Grid(runtime,data_initialization)

    #...Run the model
    # ── Timing accumulators (seconds) ────────────────────────────────────────
    PHASES = ["degradation", "uptake", "metabolism", "mortality",
              "reproduction", "output", "microbes_abundance", "repopulation"]
    t_total   = {ph: 0.0 for ph in PHASES}
    run_start = time.perf_counter()

    def _fmt(secs):
        m, s = divmod(secs, 60)
        return f"{int(m):3d}m {s:05.2f}s"

    def _print_timing(label, timings, n_days):
        phase_total = sum(timings.values())
        print(f"\n  ── Timing summary: {label} ({n_days} days) ──────────────────────────")
        for ph in PHASES:
            t = timings[ph]
            if t < 0.001:
                continue
            pct = 100 * t / phase_total if phase_total > 0 else 0
            print(f"    {ph:22s}  {_fmt(t)}  ({pct:5.1f}%)")
        print(f"    {'TOTAL (physics+output)':22s}  {_fmt(phase_total)}")
        print(f"    {'Per-day average':22s}  {phase_total/n_days*1000:.1f} ms/day", flush=True)

    for p in range(pulse):
        t_pulse = {ph: 0.0 for ph in PHASES}

        for i in range(p * cycle, (p + 1) * cycle):
            print(f"Day {i} starting", flush=True)

            t0 = time.perf_counter()
            Ecosystem.degradation(p, i)
            t_pulse["degradation"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Ecosystem.uptake(p, i)
            t_pulse["uptake"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            try:
                Ecosystem.metabolism(i)
            except Exception as e:
                print(f"CRASH in metabolism() on day {i}:", flush=True)
                traceback.print_exc()
                sys.exit(1)
            t_pulse["metabolism"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Ecosystem.mortality(i)
            t_pulse["mortality"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Ecosystem.reproduction(i)
            t_pulse["reproduction"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            if i == 0:
                Output_init.output(Ecosystem, i)
            elif i % interval == interval - 1:
                Output_init.output(Ecosystem, i)
            t_pulse["output"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Output_init.microbes_abundance(Ecosystem, i)
            t_pulse["microbes_abundance"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            if i == (p + 1) * cycle - 1:
                Ecosystem.repopulation(Output_init, i, mic_reinit)
            t_pulse["repopulation"] += time.perf_counter() - t0

        for ph in PHASES:
            t_total[ph] += t_pulse[ph]
        _print_timing(f"Pulse {p + 1}/{pulse}", t_pulse, cycle)

    # ── Full-run timing summary ───────────────────────────────────────────────
    _print_timing(f"Full run ({pulse} pulse(s))", t_total, pulse * cycle)
    print(f"    {'Wall time (incl. init)':22s}  {_fmt(time.perf_counter() - run_start)}", flush=True)

        #...export the Output_init object to the output_folder using the export() funtion in the utility module
    t0 = time.perf_counter()
    out_path = os.path.join('..', output_folder)
    os.makedirs(out_path, exist_ok=True)
    os.chdir(out_path)
    if output_format == 'csv':
        Output_init.export_to_csv(outname)
    elif output_format == 'netcdf':
        Output_init.export_to_netcdf(outname)
    elif output_format == 'zarr':
        Output_init.export_to_zarr(outname)
    else:
        warnings.warn(f"Unknown output_format '{output_format}', defaulting to NetCDF.")
        Output_init.export_to_netcdf(outname)
    print(f"\n    {'Export':22s}  {_fmt(time.perf_counter() - t0)}", flush=True)
    
main()
