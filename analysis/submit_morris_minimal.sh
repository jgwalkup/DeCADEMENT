#!/bin/bash
# Morris screening of the minimal model on HPC3 as a SLURM job array.
#
# Files needed (same layout as the repository): src/, analysis/morris_minimal.py, analysis/morris_factors_minimal.csv,
# analysis/sweep_minimal.py, this script, input_minimal_e1500_crich_km25_in/ and input/climate_data/.
# For a subset add e.g. --climates input/climate_data/climate_ambient_temp.csv --kinetics 1 to the run line below.
#
# 1. Make the design once (login node is fine, takes seconds):
#      python analysis/morris_minimal.py design morris_out --r 10 --pool 40
# 2. Submit the runs from the repository root (each array task runs every N_TASKS-th job with CPUS workers); the
#    logs/ folder must exist before sbatch (SLURM opens the --output file before the script runs):
#      mkdir -p logs
#      sbatch --array=0-$((N_TASKS-1)) analysis/submit_morris_minimal.sh morris_out
#    The array size must match N_TASKS below. Tasks are restartable: resubmit a failed index with --array=<id>.
# 3. When all tasks are done, analyse (seconds):
#      python analysis/morris_minimal.py analyze morris_out
#
# Sizing: r=10 trajectories x 27 points x 3 climates x 2 kinetics x 2 seeds = 3,240 runs of about 10 s of one core
# (20x20 grid, one 365-day pulse) = about 9-10 CPU-hours; with 32 tasks x 4 CPUs that is about 5 min of wall time per
# task plus start-up.
#SBATCH --job-name=morris_min
#SBATCH --account=allisons_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/morris_%A_%a.out

N_TASKS=32
OUT_DIR=${1:-morris_out}

module load miniconda3/25.11.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate decadementenv

mkdir -p logs
python analysis/morris_minimal.py run "$OUT_DIR" \
    --task-id "$SLURM_ARRAY_TASK_ID" --n-tasks "$N_TASKS" --procs "$SLURM_CPUS_PER_TASK"
