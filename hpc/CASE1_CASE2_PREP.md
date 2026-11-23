# BMIM Case 1 + historical Case 2: TinkerCliffs preparation

## Status, scope, and scientific warning

- Case 1 (`functional_linear`) is validated as a full 200-subject, ten-node pipeline. Its existing `R/bmim_case1_batch.R` and Slurm submission script are not replaced; one **logging-only** change fixes a bad `START:` line when `chosen_nodes` has three elements.
- Historical Case 2 (`real_or_noise`) is exposed with the *same original response generator* from `run_simulation.R`. Its **group means and SDs must be supplied from ABIDE data via a separate RDS file outside Git**. It has an independent smoke/integration check, but must still be tested on TinkerCliffs.
- **Critical:** the current Biometrics manuscript describes Case 2 using `K_gen %*% alpha` and a threshold. The original `run_simulation.R` implementation instead samples group-conditioned Gaussian responses using latent case/control labels. The repository does NOT define `K_gen` or `alpha`. We supply a separate pure function implementing the manuscript's *conditional response rule* once an authentic `K_gen` and `alpha` are provided. No kernel or alpha distribution has been invented. The legacy Case 2 results must not be presented as a numerical reproduction of the kernel-threshold equation without reconciling that scientific difference.
- The historical `real_or_noise` simulation driver did not specify a true functional interval for Case 2: `simulation_runs_new_v2.R` deliberately used `NA` domains. The new adapter records `true_domain_defined=FALSE`. **Do not calculate Case-2 domain IoU against 61:80**: those truth boundaries belong to Case 1.
- Keep `Lmax=20`, `pi_b=0.5`, strict decision threshold `P(b=1|Y)>0.5`, and the user-tested MCMC configuration of four chains, 2,500 iterations per chain, and 500 burn-in iterations unless intentionally changing the production settings.
- Inputs, RDS output, CSV exports, figures, and logs stay outside the Git repository.

## Files and exact edits

1. `R/bmim_joint_simulation.R`: minimal backward-compatible expansion. Adds `bmim_simulation_check_demo_stats()`, an optional final `demo_stats=NULL` parameter, support for `sim_setting="real_or_noise"`, and truthful domain-truth metadata. Existing Case 1 calls still pass `NA` moments as before.
2. `R/bmim_case1_batch.R`: **only the `START:` logging** is changed to print vector-valued metadata correctly. Inference, seeds, file paths, and Slurm defaults are unchanged.
3. `R/bmim_case2_batch.R`: a new restart-safe historical Case 2 array runner based on the proven Case 1 runner, requiring an external group-moments RDS file.
4. `R/bmim_case2_threshold_response.R`: a pure implementation of the manuscript's threshold response *assuming the caller explicitly provides* `K_gen`, `alpha`, and group moments. It is not wired to the historical Case 2 array, because the original code does not specify `K_gen` or `alpha`.
5. `R/checks/bmim_case2_batch_smoke.R`: CLI self-test + threshold algebra test + two tiny historical response-integration runs.
6. `hpc/bmim_case2_legacy.sbatch`: CPU-only historical Case 2 Slurm array.
7. `hpc/bmim_check_job.sh`: read-only Slurm job monitor.

## Inspect P=30 benchmark first

On TinkerCliffs login node:

```bash
squeue -u "$USER" -o '%.18i %.12P %.10T %.10M %.24R'
```

The job ID is printed by `sbatch`. When it is known:

```bash
bash "$HOME/Bayesian-Multitask-Learning/hpc/bmim_check_job.sh" JOBID
```

The results of the submitted P=30 test should be in `~/BMIM_results/case1_p30_timing/`.

## Install and test

In Windows PowerShell, from the existing project root:

```powershell
Expand-Archive -LiteralPath '.\BMIM_Case1_Case2_campaign_prep.zip' -DestinationPath '.' -Force
```

Review tracked changes before staging, since two existing files are minimally modified:

```powershell
git diff -- R/bmim_joint_simulation.R R/bmim_case1_batch.R
```

Stage **only** the listed source files:

```powershell
git add -- R/bmim_joint_simulation.R R/bmim_case1_batch.R R/bmim_case2_batch.R R/bmim_case2_threshold_response.R R/checks/bmim_case2_batch_smoke.R hpc/bmim_case2_legacy.sbatch hpc/bmim_check_job.sh hpc/CASE1_CASE2_PREP.md
```

Run `git diff --cached --check` and `git diff --cached --stat`, then commit/push. Pull with `git pull --ff-only` on TinkerCliffs **after currently running jobs finish** to keep the code version used by each scientific result unambiguous. Existing running jobs need not be cancelled.

Run CLI-only (okay on login node once R is loaded):

```bash
cd "$HOME/Bayesian-Multitask-Learning/R" && Rscript --vanilla bmim_case2_batch.R --self-test
```

On an allocated compute node, run the integrated test:

```bash
cd "$HOME/Bayesian-Multitask-Learning/R" && Rscript --vanilla -e "source('checks/bmim_case2_batch_smoke.R')"
```

Then rerun the previously passing Case 1 and joint-domain checks, because `bmim_joint_simulation.R` was extended:

```bash
cd "$HOME/Bayesian-Multitask-Learning/R" && Rscript --vanilla -e "source('checks/bmim_joint_simulation_smoke.R'); source('checks/run_all_tests_joint.R')"
```

## Prepare private group moments (historical Case 2 ONLY)

The original `simulation_runs_new_v2.R` selected `SamplingRate == 0.5` and calculated `ADOS_Total` means and SDs for PATIENT and CONTROL. It read a CSV that the user does NOT want in Git. Generate a four-number RDS summary locally, and copy **only the summary** to TinkerCliffs. The extraction requires at least two nonmissing entries per group.

PowerShell, from `C:\My Files\repos\Bayesian-Multitask-Learning` (with `$rs` bound to Rscript.exe):

```powershell
& $rs -e 'd <- read.csv("R/Demographics_ABIDE1_2_extended_F4.csv"); d <- subset(d, SamplingRate == 0.5 & !is.na(ADOS_Total)); p <- d$ADOS_Total[d$SubjectType == "PATIENT"]; c <- d$ADOS_Total[d$SubjectType == "CONTROL"]; stopifnot(length(p) >= 2, length(c) >= 2); s <- list(sim_mean_pat=mean(p), sim_sd_pat=sd(p), sim_mean_ctr=mean(c), sim_sd_ctr=sd(c)); stopifnot(all(is.finite(unlist(s))), s$sim_sd_pat > 0, s$sim_sd_ctr > 0); saveRDS(s, file.path(Sys.getenv("USERPROFILE"), "BMIM_case2_stats.rds")); print(s)'
```

On TinkerCliffs login node (the parent directory must exist before scp):

```bash
mkdir -p "$HOME/BMIM_inputs"
```

From Windows PowerShell:

```powershell
scp "$env:USERPROFILE\BMIM_case2_stats.rds" msalem@tinkercliffs1.arc.vt.edu:~/BMIM_inputs/case2_stats.rds
```

Only aggregate group summaries are transferred; no raw participant-level data are committed or copied.

## Historical Case 2 benchmark submission (NOT the manuscript kernel-threshold Case 2)

After successful integrated smoke and private summary transfer, submit one null + one signal case using a separate tag and an explicit stats path:

```bash
mkdir -p "$HOME/BMIM_results/logs" && cd "$HOME/Bayesian-Multitask-Learning" && sbatch --array=0-1%1 --export=ALL,BMIM_CASE2_STATS_RDS="$HOME/BMIM_inputs/case2_stats.rds",BMIM_N=200,BMIM_P=10,BMIM_LMAX=20,BMIM_ITER=2500,BMIM_BURNIN=500,BMIM_CHAINS=4,BMIM_PRIOR=vb,BMIM_PI_B=0.5,BMIM_REPLICATES_PER_SCENARIO=1,BMIM_RANDOMIZE_NODES=1,BMIM_BASE_SEED=31000,BMIM_RUN_TAG=case2_legacy_n200_p10_benchmark --output="$HOME/BMIM_results/logs/%x_%A_%a.out" --error="$HOME/BMIM_results/logs/%x_%A_%a.err" hpc/bmim_case2_legacy.sbatch
```

For production repeat counts, **do not silently assume 100**. Choose the required number of replicates and a concurrency cap after the P=30 timing results. For `R` replicates under each hypothesis, use array indices `0` through `2R-1` and `BMIM_REPLICATES_PER_SCENARIO=R`. The existing Case 1 batch script supports precisely that mapping too. Use different run tags and seed bases for P=10/P=30 and each response scenario.

## Safeguards

- Every task records random seed, true selected nodes, prior type, MCMC controls, and Git commit; results are saved atomically and never overwritten.
- The model and numerical sampler remain unchanged.
- `BMIM_CASE2_STATS_RDS` is mandatory; there are **no placeholder publication group means/SDs**.
- Case 2 (historical) has *no defined true time interval*, so do not report domain IoU without a new, explicit DGP specification.
- A passing smoke/integration test is not evidence of selection calibration, chain convergence, or agreement between the manuscript and historical Case-2 response mechanism.
