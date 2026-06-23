#!/bin/bash
#SBATCH --job-name=realized_traits
#SBATCH --account=allisons_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/realized_traits_%j.out
#SBATCH --error=logs/realized_traits_%j.err

# ── realized_trait_scores.py ─────────────────────────────────────────────────
# Extracts per-taxon realized trait scores (enzyme, transporter, HSP C investment)
# at substrate quartile crossing days for all communities across all scenarios.
#
# Memory notes:
#   Unlike export_summary_data.py, this script loads full taxa × days matrices
#   for each community (MicrobesSeries, Enzyme_C_Taxon, Transport_C_Taxon,
#   HSP_C_Taxon, C_Uptake_Taxon). Peak memory ≈ workers × (5 matrices × n_taxa
#   × n_days × 8 bytes). For 100 taxa × 1825 days × 5 matrices × 32 workers:
#   ~100 × 1825 × 5 × 8 × 32 ≈ 2.3 GB — 64G is conservative headroom.
#   Trim after checking seff output.
#
#   --cpus-per-task : parallel worker pool; I/O on shared HPC storage tends to
#                     saturate past ~16-32 workers, so 32 is a sensible ceiling.
#   --time          : generous bound for 9600 community runs (8 scenarios × 1200).
#                     Adjust once you have a real runtime estimate from seff.
#
# One process, many workers: keep --ntasks=1 and scale --cpus-per-task.

set -euo pipefail

PROJECT=/pub/jwalkup/DecoupledEnergyDEMENT/DecoupledEnergyDEMENTUpload
BASE=/pub/jwalkup/DecoupledEnergyDEMENT/simulations_20260606/temp_gradient
OUTCSV=/pub/jwalkup/DecoupledEnergyDEMENT/simulations_20260606_all_0615/analysis_out/realized_trait_scores.csv

cd "$PROJECT"

# Conda activation for non-interactive batch shells
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate decadementenv

mkdir -p logs
mkdir -p "$(dirname "${OUTCSV}")"

echo "Host        : $(hostname)"
echo "Start       : $(date)"
echo "CPUs alloc  : ${SLURM_CPUS_PER_TASK}"
echo "Base        : ${BASE}"
echo "Out CSV     : ${OUTCSV}"
echo

# --jobs 0 auto-detects SLURM_CPUS_PER_TASK; passed explicitly here for clarity.
python analysis/realized_trait_scores.py \
    --base    "${BASE}" \
    --out-csv "${OUTCSV}" \
    --jobs    "${SLURM_CPUS_PER_TASK}"

echo
echo "Finished    : $(date)"
