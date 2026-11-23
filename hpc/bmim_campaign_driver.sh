#!/usr/bin/env bash
# ============================================================================
# BMIM 24-hour, single-job-per-configuration campaign orchestrator.
#
# This layer changes NO model, sampler, seed, or existing result format.
# Replicates are saved by the already-validated R/bmim_case{1,2}_batch.R.
# Resumption replays task IDs 0..199; completed tasks are checked and SKIPPED
# by those R runners. Each configuration has a distinct BMIM_RUN_TAG.
#
# Normal execution: hpc/bmim_campaign.sbatch calls this script on TinkerCliffs.
# All BMIM_* settings are inherited by Slurm successors via --export=ALL.
# ============================================================================
set -Eeuo pipefail

campaign_die() { echo "CAMPAIGN ERROR: $*" >&2; exit 1; }

campaign_unsigned_integer() {
  local value="$1" name="$2" lower="$3"
  [[ "$value" =~ ^[0-9]+$ ]] || campaign_die "$name must be an integer."
  (( 10#$value >= lower )) || campaign_die "$name must be >= $lower."
}

campaign_remaining_seconds() {
  # Slurm %e reports the actual running job's expected EndTime.
  # Unknown/unparseable timestamps are fatal: never guess time left.
  local finish finish_epoch current_epoch
  finish="$(squeue --noheader --jobs "$SLURM_JOB_ID" --format='%e')" || return 1
  [[ -n "$finish" && "$finish" != "N/A" && "$finish" != "Unknown" ]] || return 1
  finish_epoch="$(date -d "$finish" +%s)" || return 1
  current_epoch="$(date +%s)" || return 1
  echo "$(( finish_epoch - current_epoch ))"
}

campaign_submit_successor() {
  [[ -z "${successor_id:-}" ]] || return 0
  (( generation < max_generations )) || {
    echo "ERROR: Maximum successor generations ($max_generations) reached." >&2
    return 1
  }
  local next_generation=$(( generation + 1 )) response
  export BMIM_CAMPAIGN_GENERATION="$next_generation"
  local -a args=(
    --parsable
    "--dependency=afterany:${SLURM_JOB_ID}"
    "--export=ALL"
    "--job-name=bmim-${case_name}-p${p}"
    "--output=${log_dir}/%x_%j.out"
    "--error=${log_dir}/%x_%j.err"
    "$repo_root/hpc/bmim_campaign.sbatch"
  )
  if ! response="$(sbatch "${args[@]}")"; then
    echo "ERROR: sbatch successor submission failed for job $SLURM_JOB_ID." >&2
    return 1
  fi
  # Slurm --parsable yields a numeric job ID, optionally followed by ;cluster.
  successor_id="${response%%;*}"
  [[ "$successor_id" =~ ^[0-9]+$ ]] || {
    echo "ERROR: Unexpected successor job ID returned by sbatch: $response" >&2
    successor_id=""
    return 1
  }
  echo "SUCCESSOR QUEUED: current=$SLURM_JOB_ID next=$successor_id dependency=afterany:$SLURM_JOB_ID generation=$next_generation"
}

campaign_on_usr1() {
  # Slurm sends this ONLY to the batch shell, not its running R child.
  # Bash processes this trap while waiting for a background Rscript.
  signal_seen=1
  echo "WALLTIME ALERT: $(date -Is) current=$SLURM_JOB_ID"
  campaign_submit_successor || echo "WARNING: Will retry successor submission after current replicate." >&2
}

campaign_on_exit() {
  local status=$?
  trap - EXIT
  # A finished campaign must not launch a redundant successor.  scancel here
  # can only target a PENDING dependent job; NEVER cancel the active job.
  if [[ "${all_complete:-0}" == "1" && -n "${successor_id:-}" ]]; then
    echo "ALL DONE: cancelling pending (not running) successor $successor_id."
    scancel "$successor_id" || echo "WARNING: Could not cancel pending successor $successor_id." >&2
  elif (( status != 0 )) && [[ -n "${successor_id:-}" ]]; then
    echo "FATAL: cancel pending successor $successor_id to avoid retry loops." >&2
    scancel "$successor_id" || true
  elif (( status == 0 )) && [[ "${all_complete:-0}" != "1" && -z "${successor_id:-}" ]]; then
    echo "FATAL: Unfinished campaign has no queued successor." >&2
    status=1
  fi
  exit "$status"
}

campaign_run_replicate() {
  local task_id="$1" pid status
  echo "REPLICATE START: $(date -Is) case=$case_name P=$p task=$task_id generation=$generation"
  BMIM_TASK_ID="$task_id" Rscript --vanilla "$runner" &
  pid=$!
  # 'wait' is interrupted by a trapped USR1. Repeat it so the R process can
  # finish normally and atomically save the result before this shell exits.
  while true; do
    if wait "$pid"; then
      status=0
      break
    else
      status=$?
      if (( status == 138 )) && (( signal_seen == 1 )); then
        signal_seen=2
        continue
      fi
      break
    fi
  done
  (( status == 0 )) || campaign_die "R runner failed for task $task_id with exit $status."
  echo "REPLICATE FINISHED: $(date -Is) case=$case_name P=$p task=$task_id"
}

campaign_main() {
  [[ -n "${SLURM_JOB_ID:-}" ]] || campaign_die "Run this via sbatch, not from the login shell."
  [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || campaign_die "This campaign uses a non-array Slurm job."
  case_name="${BMIM_CASE:-}"
  [[ "$case_name" == "case1" || "$case_name" == "case2" ]] || campaign_die "Set BMIM_CASE=case1 or case2."
  p="${BMIM_P:-}"
  [[ "$p" == "10" || "$p" == "30" ]] || campaign_die "Set BMIM_P=10 or 30."
  [[ "${BMIM_REPLICATES_PER_SCENARIO:-}" == "100" ]] || campaign_die "This production driver requires BMIM_REPLICATES_PER_SCENARIO=100."
  [[ "${BMIM_N:-}" == "200" && "${BMIM_LMAX:-}" == "20" && "${BMIM_CHAINS:-}" == "4" ]] || campaign_die "Expected N=200, LMAX=20, CHAINS=4."
  [[ "${BMIM_PRIOR:-}" == "vb" && "${BMIM_PI_B:-}" == "0.5" && "${BMIM_RANDOMIZE_NODES:-}" == "1" ]] || campaign_die "Expected vb, pi_b=0.5, randomized nodes."
  if [[ "$p" == "10" ]]; then
    [[ "${BMIM_ITER:-}" == "2500" && "${BMIM_BURNIN:-}" == "500" ]] || campaign_die "P=10 requires 2500 iterations and 500 burn-in."
    budget_seconds="${BMIM_CAMPAIGN_MIN_REMAIN_SECONDS:-900}"
  else
    [[ "${BMIM_ITER:-}" == "6000" && "${BMIM_BURNIN:-}" == "1000" ]] || campaign_die "P=30 requires 6000 iterations and 1000 burn-in."
    budget_seconds="${BMIM_CAMPAIGN_MIN_REMAIN_SECONDS:-3600}"
  fi
  campaign_unsigned_integer "$budget_seconds" BMIM_CAMPAIGN_MIN_REMAIN_SECONDS 1
  generation="${BMIM_CAMPAIGN_GENERATION:-0}"
  max_generations="${BMIM_CAMPAIGN_MAX_GENERATIONS:-30}"
  campaign_unsigned_integer "$generation" BMIM_CAMPAIGN_GENERATION 0
  campaign_unsigned_integer "$max_generations" BMIM_CAMPAIGN_MAX_GENERATIONS 1
  [[ "${BMIM_RUN_TAG:-}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || campaign_die "Set a safe, unique BMIM_RUN_TAG."
  [[ "${BMIM_BASE_SEED:-}" =~ ^[0-9]+$ ]] || campaign_die "Set BMIM_BASE_SEED explicitly."

  repo_root="${BMIM_CAMPAIGN_REPO_ROOT:-$HOME/Bayesian-Multitask-Learning}"
  [[ -f "$repo_root/hpc/bmim_campaign.sbatch" ]] || campaign_die "Missing Slurm script in $repo_root."
  if [[ "$case_name" == "case1" ]]; then
    runner="$repo_root/R/bmim_case1_batch.R"
  else
    runner="$repo_root/R/bmim_case2_batch.R"
    [[ -f "${BMIM_CASE2_STATS_RDS:-}" ]] || campaign_die "Case 2 group-moments RDS is missing."
  fi
  [[ -f "$runner" ]] || campaign_die "Missing R runner $runner."
  # Existing R runners source modules by relative paths from the R directory.
  cd "$repo_root/R"

  log_dir="${BMIM_RESULTS_DIR:-$HOME/BMIM_results}/logs"
  mkdir -p "$log_dir" "${BMIM_RESULTS_DIR:-$HOME/BMIM_results}/campaign_locks"
  local lock_file="${BMIM_RESULTS_DIR:-$HOME/BMIM_results}/campaign_locks/${BMIM_RUN_TAG}_${case_name}_p${p}.lock"
  exec 9>"$lock_file"
  flock -n 9 || campaign_die "Another campaign with this tag/case/P is running; refusing duplicate."

  # The numerical R runner performs its own CLI and numeric argument checks.
  # Do not spawn an expensive simulation until this validation succeeds.
  Rscript --vanilla "$runner" --self-test

  successor_id=""
  signal_seen=0
  all_complete=0
  trap campaign_on_usr1 USR1
  trap campaign_on_exit EXIT
  echo "CAMPAIGN START: $(date -Is) job=$SLURM_JOB_ID generation=$generation case=$case_name P=$p tasks=200 tag=$BMIM_RUN_TAG budget_seconds=$budget_seconds"
  local task_id remaining
  for (( task_id=0; task_id<200; task_id++ )); do
    remaining="$(campaign_remaining_seconds)" || campaign_die "Could not read Slurm's wall-time EndTime."
    [[ "$remaining" =~ ^-?[0-9]+$ ]] || campaign_die "Invalid time remaining: $remaining"
    if (( remaining <= budget_seconds )); then
      echo "CHECKPOINT: $remaining seconds remaining; not starting task $task_id (safe budget $budget_seconds)."
      campaign_submit_successor || campaign_die "Cannot queue required continuation."
      break
    fi
    campaign_run_replicate "$task_id"
  done

  if (( task_id == 200 )); then
    all_complete=1
    echo "CAMPAIGN COMPLETE: $(date -Is) case=$case_name P=$p all 200 tasks checked."
  else
    campaign_submit_successor || campaign_die "Unable to secure continuation."
    echo "CAMPAIGN PAUSE: $(date -Is) job=$SLURM_JOB_ID successor=$successor_id next_task=$task_id"
  fi
}

campaign_main "$@"
