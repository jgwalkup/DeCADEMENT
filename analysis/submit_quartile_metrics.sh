#!/bin/bash
#SBATCH --job-name=quartile_metrics
#SBATCH --account=allisons_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/quartile_metrics_%j.out
#SBATCH --error=logs/quartile_metrics_%j.err

# ── substrate_quartile_metrics.py  (per-community / --base mode) ─────────────────
# Each worker loads one community via compute_run, extracts its per-pulse/quartile
# crossing metrics, and returns only those small rows — so the parent never holds
# full runs and memory is bounded by the worker count, not the community count.
#   --cpus-per-task : worker pool; like the export this is I/O-bound past ~16-32
#   --mem           : ≈ workers × one run in flight + the (tiny) accumulated rows;
#                     32G is conservative headroom — trim after checking seff
#   --time          : 4000 runs (500 communities × 8 scenarios) ran the export in
#                     ~2.5 h; this loads the same files, so 6 h is a safe margin.
#                     Bump toward 8-10 h for a 1000-community (8000-run) job.

set -euo pipefail

PROJECT=/pub/jwalkup/DecoupledEnergyDEMENT/DecoupledEnergyDEMENTUpload
BASE=/pub/jwalkup/DecoupledEnergyDEMENT/simulations_20260606/temp_gradient
OUTCSV=${PROJECT}/analysis_out/substrate_quartile_metrics.csv

cd "$PROJECT"

# Conda activation for non-interactive batch shells
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate decadementenv

mkdir -p "$(dirname "$OUTCSV")"

echo "Host        : $(hostname)"
echo "Start       : $(date)"
echo "CPUs alloc  : ${SLURM_CPUS_PER_TASK}"
echo "Base        : ${BASE}"
echo "Out CSV     : ${OUTCSV}"
echo

# --jobs matched to allocated cores (the script also accepts 0 = auto-detect).
python analysis/substrate_quartile_metrics.py \
    --base    "${BASE}" \
    --out-csv "${OUTCSV}" \
    --jobs    "${SLURM_CPUS_PER_TASK}"

echo
echo "Finished    : $(date)"
