# DEMENTpy — MMRT Version

**DE**composition **M**odel of **EN**zymatic **T**raits in Python, modified to calculate
temperature-dependent reaction rates using **Macromolecular Rate Theory (MMRT)**.

A spatially explicit, trait-based model of microbial decomposition. The model simulates
a grid of microbial taxa competing for substrates via extracellular enzymes, tracking
carbon and energy flow through degradation, uptake, metabolism, mortality, and
reproduction over repeated climate pulses.

Author: **Jeth Walkup, PhD** — Department of Ecology and Evolutionary Biology,
University of California Irvine (jethwalkup@gmail.com / jwalkup@uci.edu)

## Repository layout

```
src/        Core model code (DEMENTpy engine)
  dementpy.py        Main driver / simulation loop
  initialization.py  Builds the initial ecosystem from input CSVs
  grid.py            Spatial grid; degradation, uptake, metabolism, etc.
  microbe.py         Microbial taxa and traits
  enzyme.py          Extracellular enzyme pool
  substrate.py       Substrate pool
  monomer.py         Monomer pool
  output.py          Output collection and export (CSV / NetCDF / Zarr)
  utility.py         Helper functions

input/      Model inputs (parameters, runtime config, climate data)
  runtime.txt        Simulation settings (grid size, taxa, pulses, etc.)
  parameters.csv     Model parameters
  climate.csv        Climate forcing
  ...                Additional trait / substrate / enzyme input tables

analysis/   Post-processing and figure-generation scripts
submit_multi_climate.sh   SLURM array job for HPC (UCI HPC3 / Green Planet)
```

## Requirements

Python 3 with the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

Core dependencies: numpy, pandas, scipy, xarray, matplotlib.
Optional output backends: `netCDF4` (NetCDF output) and `zarr` (Zarr output).

## Running the model

The driver reads command-line arguments and expects to be run from inside `src/`:

```
cd src
python dementpy.py <input_folder> <output_folder> <output_name> [climate_file]
```

- `input_folder` — name of the input directory (e.g. `input`)
- `output_folder` — name of the directory to write results to
- `output_name` — label for this run (digits in the name seed the RNG)
- `climate_file` — optional climate CSV (default: `climate.csv`)

Example:

```
cd src
python dementpy.py input output run001 climate.csv
```

### Configuration

Simulation settings live in `input/runtime.txt`, including grid dimensions (`x`, `y`),
number of taxa and enzymes, number of `pulse`s and time steps (`end_time`), output
`interval`, and `output_format` (`csv`, `netcdf`, or `zarr`).

## HPC / batch runs

`submit_multi_climate.sh` is a SLURM array script for running many communities across
multiple climate scenarios on a cluster. Run `bash submit_multi_climate.sh --help`-style
inspection of the header for options (community directory, climate source, partition, etc.).

## License

No license has been added yet. Without a license, default copyright applies and others
have no legal right to use, modify, or distribute this code. A license can be added later.
