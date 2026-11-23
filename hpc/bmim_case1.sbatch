#!/usr/bin/env bash
#SBATCH --job-name=bmim-case1
#SBATCH --account=iwave2
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-1%1

# A separate CPU-only Slurm allocation is created for each array task.
# These jobs survive SSH disconnection; no interactive allocation is needed.
set -Eeuo pipefail

# Executing this file must not change or terminate existing allocations.
# Output and logs live OUTSIDE Git, under $HOME/BMIM_results.
module load R/4.5.2-gfbf-2025b R-bundle-CRAN/2025.11-foss-2025b

# Keep threaded numerical libraries within this job's four allocated CPUs.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

repo_root="$HOME/Bayesian-Multitask-Learning"
if [[ ! -f "$repo_root/R/bmim_case1_batch.R" ]]; then
  echo "ERROR: batch runner missing from $repo_root/R" >&2
  exit 1
fi
cd "$repo_root/R"

echo "START: $(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-none}"
Rscript --vanilla bmim_case1_batch.R
echo "FINISHED: $(date -Is) job=${SLURM_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-none}"
