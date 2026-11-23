# BMIM joint node / functional-domain sampler (experimental, not Falcon-ready)

This archive adds **three standalone base-R modules plus one thin loader**, and **two new test scripts** to the last GP / Ising / global-Bernoulli release which passed the user's four validation scripts. The previously validated modules are retained **unchanged**; they are provided to make a reproducible local test package. Do **not** overwrite newer local source blindly.

## Mathematical model

For node `p` with `T_p` ordered functional measurements, `A_p=[s_p,e_p]` is a contiguous, **inclusive**, nonempty interval in local 1-based indices, with length `ell_p=e_p-s_p+1`, subject to:

`min_domain_length <= ell_p <= min(T_p,max_domain_length)`.

**The caller MUST set `max_domain_length`; no scientific default is chosen.** `min_domain_length` defaults to 1, and is configurable. Each node's prior is by default *uniform across admissible lengths*, and then uniform across admissible start positions conditional on its length:

`p(A_p) = [1 / (# admissible lengths)] / [T_p-ell_p+1]`.

This is NOT the old FocusGP penalty, is NOT a response-trained domain prior, and can be replaced after the experiment is validated. Node priors are the original Ising × fixed-VB empirical-Bayes weights, restricted to nonempty selections. Full functional trajectories are standardized once before MCMC. For active nodes, `D_p(A_p)` is the squared between-subject Euclidean distance averaged over retained timepoints. The GP model is

`Y | xi,A,d,g,tau2,b=1 ~ N(0,tau2 * [exp(-sum_p xi_p D_p(A_p)/d)+g*I])`.

With optional `pi_b` in (0,1), the null is `Y|b=0,sigma0^2 ~ N(0,sigma0^2 I)` with its **own** proper inverse-gamma prior. The two scale variances are integrated analytically in all node/domain/GP-model moves and redrawn afterward. Under `b=0`, the domains and indicators are auxiliary prior variables, so **interpret node-inclusion outputs conditional on `b=1` or jointly with `b=1`**, not the unconditional auxiliary xi probability.

The likelihood is conditional on fixed response-derived VB factors and on any response residualization; inference is **empirical Bayes and is not a selection-adjusted frequentist p-value**. Unlike the previous shortlist test, this sampler updates node configurations and intervals jointly without truncating their state space.

## Why node activation also proposes a domain

When `xi[p]` is switched from 0 to 1, an admissible `A[p]` is drawn from its prior. On deactivation, the same type of draw is made; the domain prior ratio cancels the forward/reverse proposal ratio exactly. This is a valid Metropolis-Hastings move on the augmented `(xi,A)` state, and prevents activation from being tied to a stale interval. Active domains also receive a mixture of independent-prior jumps and local endpoint/translation moves; the local boundary asymmetry is corrected exactly.

## Why the maximum is not a sufficient computational control

The maximum interval length reduces the candidate space and the cost of constructing each node's distance matrix. It does **not** reduce the dense GP Cholesky factorization cost `O(n_subjects^3)` by itself. By default we attempt **one domain move per sweep**, plus one node flip per node; `n_domain_updates_per_sweep` lets this be changed. All constraints remain hard throughout MCMC.

## Exact install: copy into the existing local R/ folder

**Seven files are new:**

- `bmim_domain_distances.R` — domain geometry, distance matrices, prior, and proposals
- `bmim_domain_updates.R` — active-domain and joint node/domain MH transitions
- `bmim_joint_sampler.R` — collapsed MCMC driver, global Bernoulli, scales, outputs
- `bmim_joint_domains.R` — thin source loader for the three modules
- `checks/bmim_joint_domain_smoke.R` — independent 150-state reference and detailed-balance checks
- `checks/bmim_joint_domain_gp_smoke.R` — GP/Ising/global integration smoke
- `checks/run_all_tests_joint.R` — run all previous tests plus new ones

Keep your currently validated `bmim_gp_updates.R`, `bmim_global_gp_test.R`, `bmim_global_gp_mcmc.R`, `ISGM_core.R` and the original four smoke scripts. The release includes reference copies, but none requires replacement for this extension.

### Run the independent joint toy

From `C:\My Files\repos\Bayesian-Multitask-Learning\R`:

`& $rs -e "source('checks/bmim_joint_domain_smoke.R')"`

The test enumerates all **150 discrete** `(b,xi,A)` states for 2 nodes × 5 windows each × 3 nonempty node configurations × 2 global states, with `d,g` held fixed. It checks exact prior normalization, fixed-domain identity, detailed balance of both interval and node-domain moves, and compares 4 MCMC chains to reference probabilities. It also checks the hard maximum-length constraint.

### Run the real d,g update integration smoke

`& $rs -e "source('checks/bmim_joint_domain_gp_smoke.R')"`

### Run all six scripts, including the previous four checks

`& $rs -e "source('checks/run_all_tests_joint.R')"`

**The R tests were written but could not be executed in the author's environment** (no Rscript available). The previous four scripts already passed in the user's R installation; the two new scripts still require execution on that machine. Structural delimiter checks on the new R files passed. Do not launch Falcon production runs on this unverified extension.

## Example of the joint entry point

```r
source("bmim_gp_updates.R")
source("bmim_global_gp_test.R")
source("bmim_joint_domains.R")

# Supply PRESTANDARDIZED full-time functional design and response,
# the node's full-time column lists, and fixed prior information.
# Here numbers are ILLUSTRATIVE, not recommended for the study:
fit <- bmim_joint_domain_mcmc(
  X_s = X_s,
  y_tilde = y_tilde,
  loc_slices = loc_slices,
  W_mat = W_mat,
  p_prop = p_prop,
  max_domain_length = 20L,
  n_iter = 3000L,
  burnin = 1000L,
  n_chains = 4L,
  pi_b = 0.5
)
print(fit)

# P(b=1 | Y), conditional empirical Bayes
fit$posterior_gp_included
# P(xi_p=1 | b=1,Y): node inclusion conditional on GP
fit$posterior_node_included_given_gp
# P(b=1,xi_p=1,t in A_p | Y): joint node-time GP support
fit$posterior_node_time_joint_gp
# P(xi_p=1,t in A_p | b=1,Y)
fit$posterior_node_time_given_gp
```

The helper `bmim_joint_domains_from_features(...)` reproduces the *current ISGM's* one-time `scale()` of all function columns and response, with optional nuisance-regression residualization. It requires `p_prop` explicitly to avoid silently changing the empirical-Bayes selection prior. Do not treat this as an exact drop-in to the older simulation drivers: the production drivers still need to be rewired and inspected.

## Methodological choices requiring simulation review

1. Domain prior: uniform length and uniform start, **not** the former response-tuned FocusGP criterion. If biologically informative interval-length priors are known, specify these deliberately before simulation campaigns.
2. A global indicator is optional (`pi_b=NULL` means always GP). To evaluate the Bernoulli test, specify the prior `pi_b` intentionally.
3. Conditional-EB VB probabilities `p_prop` must be constructed on the **complete trajectories** once, independent of sampled domains, and then frozen. The wrapper does not refit VB.
4. Numerical/selection calibration, scale prior sensitivity, and actual HPC throughput remain to be tested.
5. `fixed_gp` is ONLY an oracle-validation switch. Omit in full simulations.
6. The original manuscript is NOT yet updated for this joint-domain model; it must be updated after the implementation and scientific prior choices are approved.

## Test gate

Do not commit/push as production simulation-ready until the new exact-discrete validation and the complete six-test suite pass, followed by multiple-chain mixing diagnostics in a somewhat larger representative problem and a single Falcon dry run.
