#!/bin/bash
#SBATCH --job-name=export_summary
#SBATCH --account=allisons_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/export_summary_%j.out
#SBATCH --error=logs/export_summary_%j.err

# ── export_summary_data.py (streaming aggregation) ──────────────────────────────
# Memory is independent of community count, so cores are the lever, not --mem.
#   --cpus-per-task : parallel worker pool (I/O on shared storage tends to be the
#                     real ceiling past ~16-32 workers, so 32 is a sensible spot)
#   --mem           : peak ≈ (workers + 4) runs in flight + ~200 MB accumulators;
#                     48G is conservative headroom — trim it after checking seff
#   --time          : generous bound; 1200 runs (8 scenarios × 150) should finish
#                     well inside this. Adjust once you've seen one real run.
# One process, many threads: keep --ntasks=1 and scale --cpus-per-task.

set -euo pipefail

PROJECT=/pub/jwalkup/DecoupledEnergyDEMENT/DecoupledEnergyDEMENTUpload
BASE=/pub/jwalkup/DecoupledEnergyDEMENT/simulations_20260606/temp_gradient
OUTDIR=/pub/jwalkup/DecoupledEnergyDEMENT/simulations_20260606_all_0615/summary

cd "$PROJECT"

# Conda activation for non-interactive batch shells
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate decadementenv

echo "Host        : $(hostname)"
echo "Start       : $(date)"
echo "CPUs alloc  : ${SLURM_CPUS_PER_TASK}"
echo "Base        : ${BASE}"
echo "Out         : ${OUTDIR}"
echo

# --jobs 0 auto-detects SLURM_CPUS_PER_TASK; passed explicitly here for clarity.
python analysis/export_summary_data.py \
    --base    "${BASE}" \
    --out-dir "${OUTDIR}" \
    --jobs    "${SLURM_CPUS_PER_TASK}"

echo
echo "Finished    : $(date)"
