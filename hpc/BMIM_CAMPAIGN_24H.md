# BMIM 100 × 2 production campaign: 24-hour self-renewing Slurm jobs

## Scope

One initial non-array Slurm job per configuration (4 initial jobs total):

| Case | P | N | Chains | Iterations | Burn-in | Replicates H0 + H1 | Tag | Base seed |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| Case 1 (`functional_linear`) | 10 | 200 | 4 | 2500 | 500 | 100 + 100 | `production_case1_p10_100x2` | 410000 |
| Case 1 (`functional_linear`) | 30 | 200 | 4 | 6000 | 1000 | 100 + 100 | `production_case1_p30_100x2` | 420000 |
| Historical Case 2 (`real_or_noise`) | 10 | 200 | 4 | 2500 | 500 | 100 + 100 | `production_case2_p10_100x2` | 430000 |
| Historical Case 2 (`real_or_noise`) | 30 | 200 | 4 | 6000 | 1000 | 100 + 100 | `production_case2_p30_100x2` | 440000 |

All configurations use VB–Ising, randomize the three true nodes per replicate,
`Lmax=20`, `pi_b=0.5`, and strict GP classification `P(b=1|Y)>0.5`.

**Research caveat:** historical Case 2 (`real_or_noise`) is not the manuscript's
separate kernel-threshold model. The P=30 Case 2 signal benchmark had poor
mixing for `log_d` (Rhat 1.11, bulk ESS 26.4), so perform convergence checks
on the production RDS files before reporting inferential summaries.

## HPC settings and renewals

- Partition: `normal_q`, QoS: `tc_normal_base` (not the costlier 24h `tc_normal_short`).
- Job length: 24:00:00; CPUs: 4; memory: 16G; account: `iwave2`.
- No job arrays and no interactive allocations. No active allocation is
  cancelled or shortened by these scripts.
- Each job loops over task IDs 0 through 199. Tasks 0..99 are H0; 100..199 H1.
- Each result is saved atomically by the existing validated R runner; a
  completed replicate is skipped when a successor restarts the task scan.
- The orchestrator checks Slurm's actual `EndTime` before starting each new
  replicate. The default minimum safe time remaining is **15 minutes for P=10**
  and **60 minutes for P=30** (explicit conservative planning choices from
  benchmark runtimes, NOT statistical-model parameters). If less time remains,
  it does not start another expensive replicate.
- Slurm sends `SIGUSR1` to the batch shell about one hour before expiration.
  The shell immediately submits a successor via
  `sbatch --dependency=afterany:<current-job-id>`. This successor CANNOT run
  before the current job finishes, preventing overlapping writes. The in-flight
  R replicate is allowed to finish and save its RDS before the batch ends.
- If the campaign finishes after the signal, its **pending** successor is
  cancelled before it starts. No active job or interactive allocation is
  cancelled. On an R failure, any pending successor is cancelled to prevent
  infinite error loops. The maximum generation count is 30.
- The Slurm allocation ends naturally on completion or a safe checkpoint;
  no `scancel` is ever used on the running Slurm job.
- If an unexpected hard kill interrupts a replicate, that in-progress replicate
  might need recomputation. Previously saved replicates remain intact. If the
  job dies before Slurm's one-hour warning signal, automatic renewal is not
  guaranteed; inspect Slurm state/logs and resubmit manually if needed.
- Jobs must be able to invoke `sbatch` from TinkerCliffs compute nodes for
  the automatic successor feature; test this once on the live cluster.

## Install/update on TinkerCliffs

The new files are additive and do not modify the R samplers or existing job
arrays. Pull only when prior jobs using the same checkout have finished:

```bash
git -C "$HOME/Bayesian-Multitask-Learning" pull --ff-only
```

Run the **offline/mocked** smoke tests on the login node (they do NOT submit
real jobs or run R inference):

```bash
bash "$HOME/Bayesian-Multitask-Learning/hpc/checks/bmim_campaign_smoke.sh"
```

## Submit the four initial jobs (TinkerCliffs login node)

```bash
bash "$HOME/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh" case1 10
bash "$HOME/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh" case1 30
bash "$HOME/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh" case2 10
bash "$HOME/Bayesian-Multitask-Learning/hpc/bmim_campaign_submit.sh" case2 30
```

The helper checks that the real Case 2 statistics RDS exists at
`~/BMIM_inputs/case2_stats.rds`. The helper records fresh non-overlapping
base seeds and a unique run tag for each of the four configurations.
Do not submit the same configuration twice deliberately.

## Monitor

All output is outside Git, under `~/BMIM_results`:

```bash
squeue -u "$USER" -o '%.18i %.24j %.10T %.10M %.24R'
```

```bash
for tag in production_case1_p10_100x2 production_case1_p30_100x2 production_case2_p10_100x2 production_case2_p30_100x2; do d="$HOME/BMIM_results/$tag"; printf '%s: H0=%s/100 H1=%s/100\n' "$tag" "$(find "$d" -name 'truth0_rep*.rds' 2>/dev/null | wc -l)" "$(find "$d" -name 'truth1_rep*.rds' 2>/dev/null | wc -l)"; done
```

```bash
grep -hE 'CAMPAIGN START:|SUCCESSOR QUEUED:|CAMPAIGN PAUSE:|CAMPAIGN COMPLETE:|CAMPAIGN ERROR:' "$HOME"/BMIM_results/logs/bmim-*-p*_*.out | tail -n 30
```

Source files:

- `hpc/bmim_campaign.sbatch`: Slurm resources and warning signal.
- `hpc/bmim_campaign_driver.sh`: sequential task management, remaining-time
  checks, successor scheduling, exclusive lock, and signal handling.
- `hpc/bmim_campaign_submit.sh`: exact reproducible settings for each of four
  initial submissions.
- `hpc/checks/bmim_campaign_smoke.sh`: mocked end-to-end orchestration tests.

**No `.csv`, `.rds`, or other experimental data should be committed to Git.**

## VB version separation

New submissions use `BMIM_VB_SUMMARY` (`sq.corr` by default, or `mi`) and a `vb2`-specific tag. Keep any running historical campaigns and their successors on their original checkout; never pull new code into their working tree while they run. New VB batch runners refuse the old production tags and use versioned result paths. No old result is overwritten.
