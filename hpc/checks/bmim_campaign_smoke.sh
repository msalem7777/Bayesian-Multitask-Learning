#!/usr/bin/env bash
# Local/offline smoke tests for the campaign orchestration only.
# Uses mocked Slurm and R commands; NEVER submits a real job or runs inference.
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$script_dir/../.." && pwd)"
bash -n "$root/hpc/bmim_campaign.sbatch" "$root/hpc/bmim_campaign_driver.sh"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/repo/hpc" "$tmp/repo/R" "$tmp/bin" "$tmp/mock" "$tmp/results"
cp "$root/hpc/bmim_campaign.sbatch" "$tmp/repo/hpc/"
cp "$root/hpc/bmim_campaign_driver.sh" "$tmp/repo/hpc/"
touch "$tmp/repo/R/bmim_case1_batch.R" "$tmp/repo/R/bmim_case2_batch.R" "$tmp/fake_case2.rds"

cat > "$tmp/bin/Rscript" <<'MOCK'
#!/usr/bin/env bash
set -eu
for arg in "$@"; do
  if [[ "$arg" == "--self-test" ]]; then echo 'MOCK CLI: PASS'; exit 0; fi
done
echo "${BMIM_TASK_ID:-MISSING}" >> "$BMIM_MOCK_DIR/replicates"
if [[ "${BMIM_MOCK_SLEEP_FIRST:-0}" == "1" && "${BMIM_TASK_ID:-}" == "0" ]]; then sleep 1; fi
echo "MOCK REPLICATE: task=${BMIM_TASK_ID:-MISSING}"
MOCK
cat > "$tmp/bin/squeue" <<'MOCK'
#!/usr/bin/env bash
printf '%s\n' "$BMIM_MOCK_END"
MOCK
cat > "$tmp/bin/sbatch" <<'MOCK'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$BMIM_MOCK_DIR/sbatch_calls"
echo "${BMIM_CAMPAIGN_GENERATION:-UNSET}" >> "$BMIM_MOCK_DIR/generation_calls"
echo "${BMIM_CASE:-UNSET} ${BMIM_P:-UNSET} ${BMIM_ITER:-UNSET} ${BMIM_BURNIN:-UNSET} ${BMIM_REPLICATES_PER_SCENARIO:-UNSET}" >> "$BMIM_MOCK_DIR/submitted_configs"
echo '777777'
MOCK
cat > "$tmp/bin/scancel" <<'MOCK'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$BMIM_MOCK_DIR/scancel_calls"
MOCK
chmod +x "$tmp/bin/"*

export PATH="$tmp/bin:$PATH"
export BMIM_MOCK_DIR="$tmp/mock"
export BMIM_CAMPAIGN_REPO_ROOT="$tmp/repo"
export BMIM_RESULTS_DIR="$tmp/results"
export BMIM_CASE2_STATS_RDS="$tmp/fake_case2.rds"
export SLURM_JOB_ID=999900
unset SLURM_ARRAY_TASK_ID || true
export BMIM_CASE=case1 BMIM_P=10 BMIM_N=200 BMIM_LMAX=20
export BMIM_ITER=2500 BMIM_BURNIN=500 BMIM_CHAINS=4
export BMIM_PRIOR=vb BMIM_PI_B=0.5 BMIM_REPLICATES_PER_SCENARIO=100
export BMIM_RANDOMIZE_NODES=1 BMIM_BASE_SEED=123000 BMIM_RUN_TAG=smoke_case1
export BMIM_CAMPAIGN_GENERATION=0 BMIM_CAMPAIGN_MAX_GENERATIONS=30

# Test A: Enough time remaining -> run all 200 task IDs, no successor.
export BMIM_MOCK_END="$(date -d '+2 hours' +'%Y-%m-%dT%H:%M:%S')"
bash "$tmp/repo/hpc/bmim_campaign_driver.sh" > "$tmp/mock/full.log" 2>&1
[[ "$(wc -l < "$tmp/mock/replicates")" -eq 200 ]] || { echo 'FAIL: task count'; exit 1; }
[[ "$(head -n 1 "$tmp/mock/replicates")" == '0' ]] || { echo 'FAIL: first task'; exit 1; }
[[ "$(tail -n 1 "$tmp/mock/replicates")" == '199' ]] || { echo 'FAIL: last task'; exit 1; }
[[ ! -e "$tmp/mock/sbatch_calls" ]] || { echo 'FAIL: unnecessary successor'; exit 1; }
grep -q 'CAMPAIGN COMPLETE:' "$tmp/mock/full.log" || { echo 'FAIL: completion log'; exit 1; }
echo 'CAMPAIGN SMOKE, 200 sequential tasks: PASS'

# Test B: Insufficient time -> queue ONE successor, none of the R tasks.
rm -f "$tmp/mock/replicates"
export BMIM_MOCK_END="$(date -d '+5 minutes' +'%Y-%m-%dT%H:%M:%S')"
bash "$tmp/repo/hpc/bmim_campaign_driver.sh" > "$tmp/mock/renew.log" 2>&1
[[ ! -e "$tmp/mock/replicates" ]] || { echo 'FAIL: started replicate without time'; exit 1; }
[[ "$(wc -l < "$tmp/mock/sbatch_calls")" -eq 1 ]] || { echo 'FAIL: successor count'; exit 1; }
grep -q -- '--dependency=afterany:999900' "$tmp/mock/sbatch_calls" || { echo 'FAIL: successor dependency'; exit 1; }
grep -q -- '--export=ALL' "$tmp/mock/sbatch_calls" || { echo 'FAIL: env inheritance'; exit 1; }
grep -q '^1$' "$tmp/mock/generation_calls" || { echo 'FAIL: generation increment'; exit 1; }
grep -q 'CAMPAIGN PAUSE:' "$tmp/mock/renew.log" || { echo 'FAIL: pause log'; exit 1; }
echo 'CAMPAIGN SMOKE, dependent renewal: PASS'

# Test C: The same contract applies to historical Case 2 with P=30.
rm -f "$tmp/mock/sbatch_calls"
export BMIM_CASE=case2 BMIM_P=30 BMIM_ITER=6000 BMIM_BURNIN=1000 BMIM_RUN_TAG=smoke_case2
bash "$tmp/repo/hpc/bmim_campaign_driver.sh" > "$tmp/mock/case2.log" 2>&1
[[ "$(wc -l < "$tmp/mock/sbatch_calls")" -eq 1 ]] || { echo 'FAIL: case2 renewal count'; exit 1; }
grep -q -- '--job-name=bmim-case2-p30' "$tmp/mock/sbatch_calls" || { echo 'FAIL: case2 configuration'; exit 1; }
echo 'CAMPAIGN SMOKE, Case 2 P=30 routing: PASS'

# Test D: Slurm's USR1 during a running R child must enqueue a continuation,
# wait for that child, then (since all tasks finish) cancel only the pending job.
rm -f "$tmp/mock/sbatch_calls" "$tmp/mock/scancel_calls" "$tmp/mock/replicates"
export BMIM_CASE=case1 BMIM_P=10 BMIM_ITER=2500 BMIM_BURNIN=500 BMIM_RUN_TAG=smoke_signal
export BMIM_MOCK_SLEEP_FIRST=1
export BMIM_MOCK_END="$(date -d '+2 hours' +'%Y-%m-%dT%H:%M:%S')"
bash "$tmp/repo/hpc/bmim_campaign_driver.sh" > "$tmp/mock/signal.log" 2>&1 &
worker_pid=$!
for attempt in $(seq 1 200); do
  [[ -f "$tmp/mock/replicates" ]] && break
  sleep 0.01
done
[[ -f "$tmp/mock/replicates" ]] || { echo 'FAIL: signal target not launched'; exit 1; }
kill -USR1 "$worker_pid"
wait "$worker_pid" || { tail -n 40 "$tmp/mock/signal.log"; echo 'FAIL: signal killed R workflow'; exit 1; }
[[ "$(wc -l < "$tmp/mock/replicates")" -eq 200 ]] || { echo 'FAIL: signal interrupted tasks'; exit 1; }
[[ "$(wc -l < "$tmp/mock/sbatch_calls")" -eq 1 ]] || { echo 'FAIL: signal did not queue exactly one continuation'; exit 1; }
grep -q '^777777$' "$tmp/mock/scancel_calls" || { echo 'FAIL: redundant pending successor was not canceled'; exit 1; }
echo 'CAMPAIGN SMOKE, USR1 and pending-only cleanup: PASS'

# Test E: Four explicit submission commands export the exact production values.
# The mocked sbatch receives the calls; no actual Slurm jobs are submitted.
mkdir -p "$tmp/BMIM_inputs"
touch "$tmp/BMIM_inputs/case2_stats.rds"
rm -f "$tmp/mock/submitted_configs"
for case_name in case1 case2; do
  for p in 10 30; do
    HOME="$tmp" bash "$root/hpc/bmim_campaign_submit.sh" "$case_name" "$p" > "$tmp/mock/submit_${case_name}_${p}.log"
  done
done
[[ "$(wc -l < "$tmp/mock/submitted_configs")" -eq 4 ]] || { echo 'FAIL: expected 4 jobs'; exit 1; }
grep -q '^case1 10 2500 500 100$' "$tmp/mock/submitted_configs" || { echo 'FAIL: case1 P10 settings'; exit 1; }
grep -q '^case1 30 6000 1000 100$' "$tmp/mock/submitted_configs" || { echo 'FAIL: case1 P30 settings'; exit 1; }
grep -q '^case2 10 2500 500 100$' "$tmp/mock/submitted_configs" || { echo 'FAIL: case2 P10 settings'; exit 1; }
grep -q '^case2 30 6000 1000 100$' "$tmp/mock/submitted_configs" || { echo 'FAIL: case2 P30 settings'; exit 1; }
echo 'CAMPAIGN SMOKE, four initial submissions: PASS'
echo 'BMIM CAMPAIGN ORCHESTRATOR SMOKE: PASS'