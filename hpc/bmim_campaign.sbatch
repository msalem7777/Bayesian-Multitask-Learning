#!/usr/bin/env bash
#SBATCH --job-name=bmim-campaign
#SBATCH --account=iwave2
#SBATCH --partition=normal_q
#SBATCH --qos=tc_normal_base
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --signal=B:USR1@3600

# ONE Slurm job per study configuration; this job processes replicates in
# sequence.  It never cancels an active allocation or an interactive session.
# Successor jobs are queued by bmim_campaign_driver.sh only when needed.
# The --signal directive alerts the batch shell approximately one hour before
# its Slurm wall-time limit; the driver also inspects the actual EndTime.
set -Eeuo pipefail

module load R/4.5.2-gfbf-2025b R-bundle-CRAN/2025.11-foss-2025b
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

repo_root="$HOME/Bayesian-Multitask-Learning"
if [[ ! -f "$repo_root/hpc/bmim_campaign_driver.sh" ]]; then
  echo "ERROR: Missing campaign driver: $repo_root/hpc/bmim_campaign_driver.sh" >&2
  exit 1
fi

exec bash "$repo_root/hpc/bmim_campaign_driver.sh"
