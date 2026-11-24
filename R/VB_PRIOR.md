# Empirical Ising external fields (vb2)

The main joint GP likelihood, conjugate scalar effects and MCMC transitions are
unchanged. The empirical node weights are constructed once on complete
trajectories, before MCMC. The direct sampler still requires fixed `p_prop`.
The feature wrapper now accepts either explicit `p_prop` (bypasses VB) or
`p_prop=NULL` / omitted (constructs VB weights).

## Association measure

| `vb_summary` | Calculation |
|---|---|
| `auto` (default) | `mi` without Z; `residual.mi` with Z |
| `sq.corr` | Squared Pearson correlation |
| `partial.sq.corr` | Squared correlation after projecting both variables off Z |
| `mi` | Empirical mutual information of discretized measurements |
| `residual.mi` | Empirical MI after projecting both variables off Z |

MI uses the historical `infotheo::discretize` equal-frequency approach and
`mutinformation(method="emp")`, in natural-log units. We explicitly set the
number of bins to `max(2,floor(N^(1/3)))`; this documents the bin choice and
ensures at least two bins for small engineering tests. Scores are averaged
across the full time grid. Constant or numerically annihilated residual
columns/response receive zero association evidence; indices are recorded.
Scalar adjustment uses rank-aware QR with an intercept, and refuses fewer
than three residual degrees of freedom. The GP still uses the original
functional predictors and scalar-aware response likelihood, not these
residuals. Residual MI is NOT conditional MI in general. Conditional MI is
not exposed until an estimator has been agreed on and validated.

The MI defaults are candidate scientific choices, not demonstrated winners.
Finite-sample bias, binning, local signals diluted by full-domain averaging,
and interactions without marginal evidence remain important limitations.
The auxiliary Gaussian/half-normal working densities are approximations,
not the exact sampling law of either association measure.

## Auxiliary mixture and coordinate sweeps

Fit the same ClusterR two-component Gaussian mixture to the node summaries.
Use the higher component mean and variance for inclusion and the lower
component variance for a zero-centered half-normal exclusion density. Fitted
mixture proportions are not used. For identical summaries (including P=1),
no mixture separation is estimable: use zero log evidence and record status
`identical_summaries_neutral_evidence`. A failed or invalid nondegenerate
mixture stops explicitly, rather than silently substituting a different prior.

Compute Delta = log L_plus - log L_minus directly on the log scale. Initialize
m = tanh(Delta/2), then visit all P coordinates in random order each sweep.
The update is tanh(sum(W[p,]*m)+Delta[p]/2) for unit damping. The historical
one-column data-frame length bug is fixed. We retain the original default
15 sweeps and clipping to [1e-8,1-1e-8]. The objective (ELBO up to a constant)
is sum((1+m)/2 * Delta) + m' W m / 2 + Bernoulli entropy. It is evaluated
at the actual variational means after each sweep, including initialization.
The returned diagnostic history records the objective, update counts and
maximum change per sweep; fixed iterations do not guarantee convergence.

The standard `VBprop(GP=0, true_dist=0, h=0)` legacy interface delegates to
this module. Other experimental legacy branches (GP=1, true_dist=1, h=1,
and the separate VB() function) are not used by joint BMIM and have not been
re-derived or validated; their objectives are not certified by this change.
The `VBprop` scan length is corrected for the remaining branches too.

The auxiliary network is the supplied unscaled W; the GP model uses
`W * connec`, as previously. Network information enters both stages. No
frequentist calibration is claimed for the response-dependent EB prior.

## Usage (from R/)

```r
source("bmim_gp_updates.R")
source("bmim_global_gp_test.R")
source("bmim_joint_domains.R")
# X, y, W, groups and optional Z must be supplied by the caller.
fit <- bmim_joint_domains_from_features(X,W,y,groups,
  max_domain_length=20L, Z=Z, vb_summary="auto", vb_seed=123L,
  n_iter=2500L,burnin=500L,n_chains=4L,pi_b=.5)
fit$vb$association$method
fit$vb$weights
fit$vb$objective
fit$vb$updates
```

`vb_seed` (default 20261124) isolates the wrapper's random-number stream;
explicit `NULL` uses the current caller RNG. `vb_iters` defaults to 15.
Explicit `p_prop` bypasses the auxiliary calculation and preserves the
previous fixed-weight workflow. `fit$vb` is then NULL/absent. The simulation
adapter retains its replicate-specific seed and records its VB details
separately in `result$vb`. Its `prior_mode="flat"` bypasses VB entirely.

## Validation

Install `ClusterR` and `infotheo` if absent; they are required for the new
package-backed tests. From R/:

```sh
Rscript --vanilla -e "source('checks/bmim_vb_prior_smoke.R')"
Rscript --vanilla -e "source('checks/run_all_tests_joint.R')"
```

The smoke checks include exact enumeration of the variational objective,
coordinate updates, log-density stability, scalar adjustment, actual MI and
GMM, all four wrapper options, explicit weights, and reproducibility.

## HPC / old results

No-scalar batch scripts accept `BMIM_VB_SUMMARY=mi` (default) or `sq.corr`.
VB result directories now have a `vb2_mi` or `vb2_sq.corr` suffix, and the
result metadata records both algorithm version and measure. The new campaign
submitter adds the same suffix to its run tag. Old production tags are refused
by new VB batch scripts; no completed file is overwritten. Flat prior output
paths are unchanged. **Do not update the live HPC checkout while historical
jobs or their queued successors are using it.** Keep it at its existing
commit until those runs finish; use a separate checkout/worktree for new
comparisons. No job or allocation should be terminated for this amendment.
Do not pool old buggy-VB results with vb2 results. Use paired seeds and report
selection, calibration, mixing and runtime before selecting final defaults.

## Optional paired comparison

`checks/bmim_vb_comparison.R` runs a small engineering pilot with a nonlinear
signal, a scalar-confounded inactive predictor, and null scenarios. It reports
node and global Brier scores, selection metrics, interval IoU, global-indicator
switches, and elapsed times. The neutral network isolates the association
measure. Set `BMIM_VB_PILOT_REPS`, `BMIM_VB_PILOT_ITER` and (optionally)
`BMIM_VB_PILOT_OUTPUT` before sourcing. This is not the historical Case 1/2
experiment or a production convergence/calibration analysis. For paired
historical datasets, use identical batch seeds with the two explicitly named
`BMIM_VB_SUMMARY` values and separate result identities.
