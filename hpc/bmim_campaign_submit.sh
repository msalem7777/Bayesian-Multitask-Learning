#!/usr/bin/env bash
# Submit exactly ONE initial 24-hour Slurm job for a chosen BMIM configuration.
# Run on TinkerCliffs login node from anywhere:
#   bash ~/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh case1 10
#   bash ~/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh case1 30
#   bash ~/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh case2 10
#   bash ~/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh case2 30
# The running job automatically submits any needed successors itself.
set -Eeuo pipefail

if (( $# != 2 )); then
  echo 'Usage: bash hpc/bmim_campaign_submit.sh {case1|case2} {10|30}' >&2
  exit 2
fi
case_name="$1"
p="$2"
[[ "$case_name" == 'case1' || "$case_name" == 'case2' ]] || { echo 'Case must be case1 or case2.' >&2; exit 2; }
[[ "$p" == '10' || "$p" == '30' ]] || { echo 'P must be 10 or 30.' >&2; exit 2; }

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -f "$repo_root/hpc/bmim_campaign.sbatch" && -f "$repo_root/hpc/bmim_campaign_driver.sh" ]] || {
  echo 'Missing required campaign scripts. Pull the latest repository first.' >&2
  exit 1
}

# Scientific settings validated in the four full-size benchmarks.
if [[ "$p" == '10' ]]; then
  iterations=2500
  burnin=500
else
  iterations=6000
  burnin=1000
fi
if [[ "$case_name" == 'case1' && "$p" == '10' ]]; then base_seed=410000; fi
if [[ "$case_name" == 'case1' && "$p" == '30' ]]; then base_seed=420000; fi
if [[ "$case_name" == 'case2' && "$p" == '10' ]]; then base_seed=430000; fi
if [[ "$case_name" == 'case2' && "$p" == '30' ]]; then base_seed=440000; fi

# New VB algorithms must never share historical campaign result identities.
vb_summary="${BMIM_VB_SUMMARY:-mi}"
[[ "$vb_summary" == 'mi' || "$vb_summary" == 'sq.corr' ]] || { echo 'BMIM_VB_SUMMARY must be mi or sq.corr.' >&2; exit 2; }
tag="production_${case_name}_p${p}_100x2_vb2_${vb_summary}"
results_dir="$HOME/BMIM_results"
mkdir -p "$results_dir/logs"

if [[ "$case_name" == 'case2' ]]; then
  stats_path="$HOME/BMIM_inputs/case2_stats.rds"
  [[ -r "$stats_path" ]] || {
    echo "ERROR: Missing Case 2 statistics: $stats_path" >&2
    exit 1
  }
  export BMIM_CASE2_STATS_RDS="$stats_path"
fi

# Slurm gives preexisting environment values precedence over assignments in
# --export=ALL,NAME=value. Explicitly export the desired scientific controls
# before sbatch to prevent accidentally inheriting stale BMIM_* settings.
export BMIM_CASE="$case_name" BMIM_N=200 BMIM_P="$p" BMIM_LMAX=20
export BMIM_ITER="$iterations" BMIM_BURNIN="$burnin" BMIM_CHAINS=4
export BMIM_VB_SUMMARY="$vb_summary"
export BMIM_PRIOR=vb BMIM_PI_B=0.5 BMIM_REPLICATES_PER_SCENARIO=100
export BMIM_RANDOMIZE_NODES=1 BMIM_BASE_SEED="$base_seed" BMIM_RUN_TAG="$tag"
export BMIM_CAMPAIGN_GENERATION=0 BMIM_CAMPAIGN_MAX_GENERATIONS=30
export BMIM_RESULTS_DIR="$results_dir"
if [[ "$p" == '10' ]]; then
  export BMIM_CAMPAIGN_MIN_REMAIN_SECONDS=900
else
  export BMIM_CAMPAIGN_MIN_REMAIN_SECONDS=3600
fi

# Do not use --array: this is one job processing 200 replicates sequentially.
# --export=ALL preserves the exact BMIM_* variables in successor jobs.
job_id="$(sbatch --parsable \
  --job-name="bmim-${case_name}-p${p}" \
  --output="$results_dir/logs/%x_%j.out" \
  --error="$results_dir/logs/%x_%j.err" \
  --export=ALL \
  "$repo_root/hpc/bmim_campaign.sbatch")"

echo "SUBMITTED: job=${job_id} case=${case_name} P=${p} N=200 null=100 signal=100 time=24:00:00 tag=${tag}"
echo "LOGS: $results_dir/logs/bmim-${case_name}-p${p}_<jobid>.out"
echo 'SUCCESSORS: automatic, with afterany dependency and per-replicate checkpointing.'
