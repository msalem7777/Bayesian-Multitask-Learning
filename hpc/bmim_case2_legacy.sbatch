#!/usr/bin/env bash
#SBATCH --job-name=bmim-case2-legacy
#SBATCH --account=iwave2
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-1%1

# The historical R Case 2 response mixture is distinct from the manuscript's
# unimplemented K_gen %*% alpha threshold model. Do not silently conflate them.
# No dataset or result files are ever placed in the Git repository.
set -Eeuo pipefail
module load R/4.5.2-gfbf-2025b R-bundle-CRAN/2025.11-foss-2025b
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

repo_root="$HOME/Bayesian-Multitask-Learning"
cd "$repo_root/R"
if [[ -z "${BMIM_CASE2_STATS_RDS:-}" || ! -f "${BMIM_CASE2_STATS_RDS}" ]]; then
  echo "ERROR: BMIM_CASE2_STATS_RDS must name a readable external group-moments RDS" >&2
  exit 1
fi

echo "START: $(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-none}"
Rscript --vanilla bmim_case2_batch.R
echo "FINISHED: $(date -Is) job=${SLURM_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-none}"
