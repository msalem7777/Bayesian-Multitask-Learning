#!/usr/bin/env bash
# Read-only dashboard for a Slurm array. Usage: bash hpc/bmim_check_job.sh JOBID
set -Eeuo pipefail
job_id="${1:-}"
if [[ ! "$job_id" =~ ^[0-9]+$ ]]; then
  echo 'Usage: bash hpc/bmim_check_job.sh <numeric-array-job-id>' >&2
  exit 2
fi
printf '\n=== Pending/running ===\n'
squeue -j "$job_id" -o '%.18i %.12P %.10T %.10M %.24R' || true
printf '\n=== Accounting (may lag) ===\n'
sacct -j "$job_id" --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode || true
printf '\n=== Recent task log lines ===\n'
shopt -s nullglob
logs=("$HOME"/BMIM_results/logs/*_"${job_id}"_*.out)
if (( ${#logs[@]} == 0 )); then
  echo 'No task logs found yet (pending tasks may not have logs).'
else
  for logfile in "${logs[@]}"; do
    echo "--- ${logfile} ---"
    tail -n 6 "$logfile"
  done
fi
