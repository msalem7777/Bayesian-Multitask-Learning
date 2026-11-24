# Scalar regression effects in joint BMIM

## Model and defaults

When a scalar design `Z` is supplied, the original response is modeled jointly
with the functional effect, not residualized by OLS. Let `C = K + g I` and
let `B` be a positive definite coefficient prior covariance multiplier.

- GP alternative: `y | beta,tau2 ~ N(Z beta, tau2 C)` and
  `beta | tau2,g ~ N(mu, tau2 g B)`.
- Noise-only null: `y | beta,sigma0^2 ~ N(Z beta, sigma0^2 I)` and
  `beta | sigma0^2 ~ N(mu, sigma0^2 B)`.
- The model-specific inverse-gamma scale priors, GP priors, Ising/VB weights,
  and domain priors otherwise remain unchanged. `mu = 0` and `B = 100 I` are
  configurable engineering defaults, not calibrated scientific choices.
  In particular, the prior coefficient SD is ten times the noise SD.
- The same beta prior conditional on observation-noise variance is used in
  both models. The alternative and null retain distinct variance scales.

### Default: `beta_mode = "integrated"`

Integrate beta exactly before integrating the active variance. Alternative
covariance multiplier is `K + g*(I + Z B Z^T)` and the null multiplier is
`I + Z B Z^T`, with mean `Z mu` in both cases. `Z B Z^T` is computed once.
The usual inverse-gamma integration uses response dimension N. All proposals
use this same likelihood. `save_beta=TRUE` optionally generates Gaussian
conditional draws after each sweep; those draws do not enter subsequent
transitions (but consume RNG draws). The default saves no beta samples.

### Optional: `beta_mode = "sampled"`

Keep a single shared beta in the chain. Hold beta fixed while updating
`d,g,xi,A,b`, integrating the active variance in every move. For the GP
alternative, use the joint density of y and beta: determinant term
`log|C| + log|B| + S log(g)` and quadratic term
`(y-Z beta)^T C^-1 (y-Z beta) + (beta-mu)^T B^-1 (beta-mu)/g`.
The integrated variance shape is `a + (N+S)/2`, not `a + N/2`.
Under the null use `C=I`, `g=1`, and the null IG prior. Including the beta
prior density in the global-indicator and nugget updates is essential.
After the collapsed moves draw the active variance conditional on beta, then
sample beta from its exact multivariate Gaussian conditional. The inactive
variance is refreshed from its prior. Both modes thus target the same
conditional empirical-Bayes posterior for common priors and data; their
mixing and finite-chain estimates can differ. Sampled mode always stores beta.

## Usage

From the repository `R/` directory, source the existing loader sequence:

```r
source("bmim_gp_updates.R")
source("bmim_global_gp_test.R")
source("bmim_joint_domains.R")

# Use your actual functional data, unresidualized response, scalar matrix,
# fixed inclusion weights and interval limit. MCMC settings are examples.
fit <- bmim_joint_domains_from_features(
  path_features = X, W_mat = W, y = y, grps_vec = groups,
  p_prop = p_prop, max_domain_length = 20L, Z = scalar_covariates,
  n_iter = 3000L, burnin = 1000L, n_chains = 4L, pi_b = 0.5,
  beta_mode = "integrated", save_beta = TRUE
)
# Same API for explicit sampling: replace beta_mode with "sampled".
fit$beta_mode
fit$beta_prior_cov
fit$response_preprocessing
fit$chains[[1]]$beta # NULL unless requested, or always present in sampled mode
```

`beta_prior_mean` is a numeric vector and `beta_prior_cov` is an S-by-S
positive definite matrix in the coefficient coordinates actually used by
the sampler. At the wrapper level these include an intercept first, followed
by centered/scaled scalar covariates. The wrapper rejects constant scalar
columns; encode categorical variables numerically before calling it. Use a
zero-column N-row numeric matrix for `Z` to fit an intercept only. At the
lower-level `bmim_joint_domain_mcmc(..., Z=...)`, the supplied design is used
exactly, with no added intercept or standardization. The historical argument
`y_tilde` must receive the unresidualized response when `Z` is supplied.

**Response scaling:** with scalar covariates, the wrapper uses original response
units by default, or divides by a positive fixed `response_scale` specified
independently of the observed response. It does not center or estimate a
response SD. Inverse-gamma rates are in the chosen response-variance units;
review them for the scientific application. The wrapper records its scaling
and scalar-column centering/scaling for interpreting beta draws. With no
scalar design and no fixed response scale the wrapper preserves the historical
standardized-response likelihood and existing no-covariate campaign. A numeric
`response_scale` also opts out of historical response standardization when
Z is NULL; no intercept is automatically added in that case. The direct sampler
performs no response scaling. Do not compare evidence between fits with
inconsistent preprocessing/prior choices.

The caller supplies fixed `p_prop`, exactly as before. No VB weights are
recomputed internally. Whether these weights should be constructed using
scalar-adjusted associations requires a separate scientific decision. The
resulting global probability is conditional empirical Bayes, not a calibrated
frequentist p-value. Existing fixed-domain and historical ABIDE scripts are
not silently converted; call the joint wrapper for scalar-aware analysis.

## Validation gate

Run the standalone scalar check first from `R/`:

```sh
Rscript --vanilla -e "source('checks/bmim_scalar_effects_smoke.R')"
```

It checks independent dense Gaussian/inverse-gamma evidence, joint density
constants with unequal null/alternative variance priors, Gaussian conditional
moments, finite-state enumeration with both coefficient modes, variable d/g
with slice and random-walk proposals, selection-only and global models,
wrapper preprocessing, and no-design compatibility. The finite-state checks
are stochastic, not a production mixing or calibration guarantee.

Then run the full regression suite:

```sh
Rscript --vanilla -e "source('checks/run_all_tests_joint.R')"
```

No R executable is available in the artifact-generation environment and the
attempt to access the Debian package index failed because DNS was unavailable.
These R tests are therefore written but not executed. Do not treat the scalar
update as production-validated until these tests pass locally. No simulation or
ABIDE results have been generated. The default coefficient prior and variance
scaling require scientific review before manuscript reporting.
