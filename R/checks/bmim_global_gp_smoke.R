# ============================================================================
# Smoke test for the CONDITIONAL BMIM global indicator comparison.
# From the R directory:
#   & $rs -e "source('checks/bmim_global_gp_smoke.R')"
# No external packages. This file does not run production ISGM.
# ============================================================================
script_files <- Filter(Negate(is.null),
                       lapply(sys.frames(), function(frame) frame$ofile))
if (!length(script_files)) stop("Source this test from a real file.")
r_dir <- dirname(dirname(normalizePath(tail(script_files, 1L)[[1L]])))
source(file.path(r_dir, "bmim_gp_updates.R"))
source(file.path(r_dir, "bmim_global_gp_test.R"))

# Load the global Bernoulli sampler before checking its
# conditional noise-only Gaussian likelihood.
source(file.path(r_dir, "bmim_global_gp_mcmc.R"))

# 1. Small deterministic model with three nodes and seven allowed states.
set.seed(20261122)
n <- 9L
z <- seq(-1, 1, length.out = n)
X <- scale(cbind(z, z^2, sin(3*z), cos(z), z^3, sin(5*z)))
y <- as.numeric(scale(0.6*X[,1] + 0.3*X[,3] + sin(z*2)))
D <- bmim_gp_precompute_distances(X, list(1:2, 3:4, 5:6))
prior <- bmim_gp_default_prior(D)
W <- matrix(c(0, .2, .1, .2, 0, .3, .1, .3, 0), 3, 3)
r <- c(.67, .35, .55)
config <- as.matrix(expand.grid(rep(list(0:1), 3)))
config <- config[rowSums(config) > 0, , drop = FALSE]
storage.mode(config) <- "double"
stopifnot(nrow(config) == 7L)

# Construct an illustrative empirical selection trace with known frequencies;
# this intentionally bypasses MCMC to independently test shortlist logic.
freq <- c(35L, 25L, 15L, 10L, 7L, 5L, 3L)
trace <- config[rep(seq_len(7L), freq), , drop = FALSE]
s <- bmim_shortlist(trace, M = 3L, threshold = .05, mass_target = .90)
stopifnot(nrow(s$configurations) == 5L)
stopifnot(abs(s$covered_mass - .92) < 1e-12)
stopifnot(abs(s$omitted_mass - .08) < 1e-12)
cat("SHORTLIST RULE AND COVERAGE: PASS\n")

# 2. Null sigma0^2 has an independent, identifiable IG prior.
# The null evidence has no dependence on the GP nugget, GP scale, or d.
a0 <- prior$null_sigma_shape
rate0 <- prior$null_sigma_rate
manual <- lgamma(a0+n/2) - lgamma(a0) + a0*log(rate0) -
  n*log(2*pi)/2 - (a0+n/2)*log(rate0 + sum(y*y)/2)
stopifnot(isTRUE(all.equal(
  manual, bmim_global_null_logevidence(y, prior),
  tolerance = 1e-10)))
# Conditional Gaussian log-likelihood uses only sigma0^2.
sigma2_value <- 0.8
stopifnot(isTRUE(all.equal(
  bmim_global_null_loglik(y, sigma2_value),
  -0.5*(n*log(sigma2_value)+sum(y*y)/sigma2_value),
  tolerance = 1e-12)))
# Changing the null prior changes only the null evidence, not GP likelihood.
prior_null_shift <- prior
prior_null_shift$null_sigma_rate <- 2*prior$null_sigma_rate
stopifnot(abs(bmim_global_null_logevidence(y, prior_null_shift) -
              manual) > 1e-5)
cat("INDEPENDENT NULL-VARIANCE MARGINAL LIKELIHOOD: PASS\n")

# 3. Compare shortlisting and full seven-state evidence on THE SAME MC draws.
# This checks algebra independent of random integration convergence.
set.seed(444)
draws <- bmim_global_draw_prior(prior, 250L)
args <- list(y_tilde=y, D_by_loc=D, gp_prior=prior,
             W_mat=W, p_prop=r, pi_b=.5, prior_draws=draws)
full <- do.call(bmim_global_gp_test,
                c(list(selection_samples=trace, M=7L,
                       threshold=0, mass_target=1), args))
short <- do.call(bmim_global_gp_test,
                 c(list(selection_samples=trace, M=3L,
                        threshold=.05, mass_target=.90), args))
stopifnot(nrow(full$configurations) == 7L)
stopifnot(nrow(short$configurations) == 5L)
stopifnot(is.finite(full$log_bayes_factor_10_truncated))
stopifnot(is.finite(short$log_bayes_factor_10_truncated))
stopifnot(isTRUE(all.equal(full$log_evidence_null,
                           short$log_evidence_null, tolerance = 1e-12)))

# Independent reconstruct: weighted sum over ALL alternatives.
w <- bmim_global_log_model_prior(full$configurations, W, r)
independent_log_alt <- bmim_logsumexp(full$log_evidence_each_config + w)
stopifnot(isTRUE(all.equal(independent_log_alt,
                           full$log_evidence_alternative_truncated,
                           tolerance=1e-12)))
stopifnot(isTRUE(all.equal(
  full$posterior_gp_included_truncated,
  plogis(full$log_bayes_factor_10_truncated), tolerance=1e-12)))
cat("FULL-SEVEN AND TRUNCATED BAYES ALGEBRA: PASS\n")

# 4. Confirm all original configuration prior weights sum to 1 on shortlist.
stopifnot(abs(sum(short$conditional_model_prior)-1) < 1e-12)
stopifnot(abs(sum(full$conditional_model_prior)-1) < 1e-12)
stopifnot(short$posterior_gp_included_truncated >= 0 &&
          short$posterior_gp_included_truncated <= 1)
stopifnot(abs(sum(full$conditional_model_posterior_given_gp)-1) < 1e-12)
set.seed(334)
indicator_draws <- bmim_global_draw_indicator(short, n=200L)
stopifnot(nrow(indicator_draws) == 200L,
          all(indicator_draws$b %in% c(0L, 1L)),
          all(is.na(indicator_draws$configuration[indicator_draws$b == 0L])),
          all(!is.na(indicator_draws$configuration[indicator_draws$b == 1L])))
cat("BINARY MODEL PRIOR AND POSTERIOR: PASS\n")

# 5. Test the convenience path's burn-in handling without fitting ISGM.
details <- list(samples = rbind(config[1,,drop=FALSE],trace),
                y_tilde=y, D_by_loc=D, gp_prior=prior, W_mat=W, p_prop=r)
via_details <- bmim_global_gp_from_isgm(details, burnin=1L,
                                      M=7L, mass_target=1, pi_b=.5,
                                      prior_draws=draws)
stopifnot(isTRUE(all.equal(via_details$log_bayes_factor_10_truncated,
                           full$log_bayes_factor_10_truncated,
                           tolerance=1e-12)))
cat("OPTIONAL ISGM DETAILS INTERFACE: PASS\n")
# 6. Actual ISGM integration (short run, no VB fit required).
source(file.path(r_dir, "ISGM_core.R"))
set.seed(777)
details_live <- ISGM(path_features = X, W_mat = W, y = y,
                     curr_xi = c(1, 1, 0),
                     grps_vec = list(list(1:2), list(1:2), list(1:2)),
                     N = 8L, VBprop = 0, return_details = TRUE)
stopifnot(is.list(details_live), nrow(details_live$samples) == 8L,
          ncol(details_live$samples) == 3L,
          all(rowSums(details_live$samples) >= 1L),
          length(details_live$p_prop) == 3L)
# Confirm that the non-opt-in legacy call remains a simple N x P matrix.
set.seed(777)
legacy_live <- ISGM(path_features = X, W_mat = W, y = y,
                    curr_xi = c(1, 1, 0),
                    grps_vec = list(list(1:2), list(1:2), list(1:2)),
                    N = 8L, VBprop = 0)
stopifnot(is.matrix(legacy_live), nrow(legacy_live) == 8L,
          isTRUE(all.equal(unname(details_live$samples),
                           unname(legacy_live), tolerance=0)))
cat("ISGM MATRIX COMPATIBILITY AND DETAILS INTEGRATION: PASS\n")

# 7. True Bernoulli/GP MCMC smoke with distinct initial b states.
result_mcmc <- bmim_global_gp_mcmc_from_isgm(
  details_live, burnin_selection=1L, M=2L, pi_b=.5,
  n_iter=40L, burnin=10L, n_chains=2L,
  threshold=.05, mass_target=.90, seed=123L)
stopifnot(inherits(result_mcmc, "bmim_global_gp_mcmc_result"),
          length(result_mcmc$chains)==2L,
          all(vapply(result_mcmc$chains, nrow, integer(1))==40L),
          nrow(result_mcmc$samples_after_burnin)==60L,
          all(result_mcmc$samples_after_burnin[, "b"] %in% c(0,1)),
          all(result_mcmc$samples_after_burnin[, "tau2"] > 0),
          all(result_mcmc$samples_after_burnin[, "sigma2"] > 0),
          all(is.finite(result_mcmc$samples_after_burnin)),
          abs(sum(result_mcmc$joint_posterior_gp_configurations)-
              result_mcmc$posterior_gp_included_truncated)<1e-12)
cat("ACTUAL JOINT BERNOULLI/GP MCMC SMOKE: PASS\n")

cat("BMIM CONDITIONAL GLOBAL-GP TESTS: PASS\n")
cat("Note: this is a conditional data-adaptive EB comparison, not an exact\n")
cat("      post-selection p-value or integrated full-BMIM Bayes factor.\n")
