# BMIM Case 1 Slurm batch execution

This adds a CPU-only job-array runner for the **functional-linear (Case 1)**
scenario already exercised by `bmim_joint_simulation.R`. It does not modify
GP/Ising/VB/joint-domain inference code, `run_simulation.R`, or any validated
scripts. It does **not** yet implement Case 2 or the old competing methods.

## Defaults are a benchmark, not publication settings

The two default array tasks are one noise-only dataset (task 0) and one
signal-present dataset (task 1), with **N=200, P=10, Lmax=20, 1000 iterations,
250 burn-in, 4 chains, VB, and pi_b=0.5**. The decision is strictly
`P(b=1 | Y) > 0.5`.

Active nodes are `2,5,8` by default, matching the reference script's Case 1
configuration. To randomize active-node identities per replicate, set
`BMIM_RANDOMIZE_NODES=1` explicitly; this is necessary if we want to match the
original paper's claim that they vary from replicate to replicate.

## On Windows: install, stage, and push

Extract the provided ZIP in your local repository root. The new files are:

- `R/bmim_case1_batch.R`
- `R/checks/bmim_case1_batch_smoke.R`
- `hpc/bmim_case1.sbatch`
- `hpc/README.md`

Only stage those four files. Keep CSV/PNG/RDS results out of GitHub.

## On TinkerCliffs: pull and smoke-test on the current compute allocation

```bash
cd "$HOME/Bayesian-Multitask-Learning" && git pull --ff-only && cd R && Rscript --vanilla -e "source('checks/bmim_case1_batch_smoke.R')"
```

## Submit an SSH-independent two-task benchmark

From **any TinkerCliffs login shell**, run:

```bash
mkdir -p "$HOME/BMIM_results/logs" && cd "$HOME/Bayesian-Multitask-Learning" && sbatch --output="$HOME/BMIM_results/logs/%x_%A_%a.out" --error="$HOME/BMIM_results/logs/%x_%A_%a.err" hpc/bmim_case1.sbatch
```

Slurm allocates new CPU resources for each array task. It does **not** use,
modify, or terminate any running interactive allocation or Falcon GPU job.
`--array=0-1%1` is the script's safe default: two tasks, sequential execution.

The output files are under
`$HOME/BMIM_results/pilot/case1_N200_P10_L20/iter1000_burn250_ch4_vb_pi0.5_random0/`.
The output for each task is an `.rds` containing full fitted chains, truth,
prior weights, network, configuration, source commit, and runtime. Results are
written to a temporary sibling and renamed after successful completion.
Rerunning the same task skips a completed RDS instead of overwriting it.

## To launch repeated simulations AFTER approving all scientific settings

`BMIM_REPLICATES_PER_SCENARIO=R` means R null and R alternative replicates.
The task IDs must be exactly `0,...,2*R-1`; task IDs less than R are null,
remaining task IDs are alternative. Choose a new run tag to avoid accidental
reuse of a pilot result.

For example, **this is an illustration, not an approved production request**:

```bash
mkdir -p "$HOME/BMIM_results/logs" && cd "$HOME/Bayesian-Multitask-Learning" && sbatch --array=0-19%2 --export=ALL,BMIM_REPLICATES_PER_SCENARIO=10,BMIM_RUN_TAG=case1_production_review,BMIM_N=200,BMIM_P=10,BMIM_LMAX=20,BMIM_ITER=4000,BMIM_BURNIN=1000,BMIM_CHAINS=4,BMIM_PRIOR=vb,BMIM_PI_B=0.5,BMIM_RANDOMIZE_NODES=1 --output="$HOME/BMIM_results/logs/%x_%A_%a.out" --error="$HOME/BMIM_results/logs/%x_%A_%a.err" hpc/bmim_case1.sbatch
```

The numeric settings above are illustrative only. Do not submit the production
example without deciding the Monte Carlo replication count, MCMC length,
resource limits, and whether to randomize node identities.

## Environment variables

| Variable | Pilot default | Meaning |
|---|---|---|
| `BMIM_N` | `200` | Subjects |
| `BMIM_P` | `10` | Functional nodes |
| `BMIM_LMAX` | `20` | Hard maximum contiguous interval length |
| `BMIM_ITER` | `1000` | Iterations in each chain |
| `BMIM_BURNIN` | `250` | Initial discarded iterations |
| `BMIM_CHAINS` | `4` | Independent chains |
| `BMIM_REPLICATES_PER_SCENARIO` | `1` | Number of null AND number of signal replicates |
| `BMIM_PRIOR` | `vb` | `vb` or `flat` |
| `BMIM_PI_B` | `0.5` | Prior probability the GP is present |
| `BMIM_BASE_SEED` | `19700` | Seed of task 0; task ID is added |
| `BMIM_RUN_TAG` | `pilot` | Results directory name; set uniquely per study |
| `BMIM_CHOSEN_NODES` | `2,5,8` | Fixed true signal nodes |
| `BMIM_RANDOMIZE_NODES` | `0` | Set `1` to vary three nodes by replicate |
| `BMIM_RESULTS_DIR` | `$HOME/BMIM_results` | External result storage |

## Job monitoring

```bash
squeue -u "$USER" -o '%.18i %.12P %.10T %.10M %.24R'
```

The submitted command prints its Slurm job ID. The job persists after SSH
ends. A Slurm batch job can still fail if time/memory limits are exceeded.
When a task fails, its incomplete temporary `.rds` is not reported as a
completed result.

## Other checks

From `R/`:

```bash
Rscript --vanilla -e "source('checks/bmim_case1_batch_smoke.R')"
```

For the previously passing seven checks:

```bash
Rscript --vanilla -e "source('checks/run_all_tests_joint.R'); source('checks/bmim_joint_simulation_smoke.R'); source('checks/bmim_simulation_covariance_smoke.R')"
```

The batch runner checks for the covariance repair, requires ClusterR in VB
mode, and records the strict global decision `P(GP) > 0.5`.
