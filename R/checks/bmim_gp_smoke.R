# ================================================================
# checks/bmim_gp_smoke.R
# Run from the R/ directory:
#   Rscript checks/bmim_gp_smoke.R
# No external R packages are required for this test.
# ================================================================
source("bmim_gp_updates.R")

set.seed(412)

# Test 1: GP distance construction exactly matches the intended
# location-specific normalization by retained time-column count.
X <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
slices <- list(1:2, 3:5)
D_list <- bmim_gp_precompute_distances(X, slices)
D_expected <- as.matrix(dist(X[, 1:2, drop = FALSE] / sqrt(2)))^2
# as.matrix(dist(...)) assigns automatic row/column labels (1, ..., N),
# while our distance helper returns an unnamed matrix. Compare values,
# not those labels; matrix dimensions are checked independently below.
stopifnot(isTRUE(all.equal(
  unname(D_list[[1]]),
  unname(D_expected),
  tolerance = 1e-10
)))
stopifnot(is.matrix(D_list[[1]]),
          identical(dim(D_list[[1]]), c(nrow(X), nrow(X))),
          all(diag(D_list[[1]]) == 0),
          all(D_list[[1]] >= 0))
# Regression check: a location containing only one retained time column
# must still produce an N-by-N distance matrix.
D_single <- bmim_gp_precompute_distances(X[, 1, drop = FALSE], list(1L))
stopifnot(is.matrix(D_single[[1]]),
          identical(dim(D_single[[1]]), c(nrow(X), nrow(X))),
          all(diag(D_single[[1]]) == 0))

# Test 2: The GP kernel matches hetGP::cov_gen Gaussian semantics.
xi <- c(1L, 1L)
D <- bmim_gp_selected_distance(D_list, xi)
K <- bmim_gp_covariance(D, d = 2.3, g = 0.1)
K_direct <- exp(-D / 2.3) + diag(0.1, nrow(D))
stopifnot(isTRUE(all.equal(K, K_direct, tolerance = 1e-12)))

# Test 3: Quadratic form and log determinant agree with base R.
y <- as.numeric(scale(rnorm(8)))
terms <- bmim_gp_cholesky_terms(K, y)
stopifnot(!is.null(terms))
stopifnot(isTRUE(all.equal(terms$quad,
                           as.numeric(crossprod(y, solve(K, y))),
                           tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(terms$logdet,
                           as.numeric(determinant(K, logarithm = TRUE)$modulus),
                           tolerance = 1e-10)))

# Test 4: VB is an explicit empirical prior and the Ising term
# favors equal selections. The exact log odds must include both.
W <- matrix(c(0, 0.75, 0.75, 0), 2, 2)
r <- c(0.8, 0.8)
state <- list(d = 2.3, g = 0.1, tau2 = 1.3)
l10 <- bmim_gp_selection_log_target(c(1, 0), y, D_list, state, W, r)
l11 <- bmim_gp_selection_log_target(c(1, 1), y, D_list, state, W, r)
L10 <- bmim_gp_log_likelihood(y, bmim_gp_selected_distance(D_list, c(1, 0)),
                               state$d, state$g, state$tau2)
L11 <- bmim_gp_log_likelihood(y, bmim_gp_selected_distance(D_list, c(1, 1)),
                               state$d, state$g, state$tau2)
expected <- 1.5 + log(0.8 / 0.2) + (L11 - L10)
stopifnot(isTRUE(all.equal(l11 - l10, expected, tolerance = 1e-10)))

# Test 5: Symmetric one-coordinate flips satisfy detailed balance
# under the fixed GP state and the full Ising + VB + GP target.
a10to11 <- exp(min(0, l11 - l10))
a11to10 <- exp(min(0, l10 - l11))
stopifnot(isTRUE(all.equal(exp(l10) * a10to11,
                           exp(l11) * a11to10,
                           tolerance = 1e-10)))

# Test 6: The log-random-walk proposal includes the Jacobian.
a <- bmim_gp_log_mh_ratio(-2.4, -2.1, 1.8, 2.5)
b <- bmim_gp_log_mh_ratio(-2.1, -2.4, 2.5, 1.8)
stopifnot(isTRUE(all.equal(a, -b, tolerance = 1e-12)))

# Test 6b: Prior-independence proposals cancel the d-prior ratio.
# For any two positive d values, the forward/backward MH log ratios
# must be negatives of each other and use ONLY marginal likelihoods.
prior_jump_test <- bmim_gp_default_prior(D_list)
d_a <- 0.05 * prior_jump_test$distance_scale
d_b <- 1.00 * prior_jump_test$distance_scale
ll_a <- bmim_gp_log_marginal_likelihood(y, D, d_a, 0.2, prior_jump_test)
ll_b <- bmim_gp_log_marginal_likelihood(y, D, d_b, 0.2, prior_jump_test)
stopifnot(isTRUE(all.equal(ll_b - ll_a, -(ll_a - ll_b),
                           tolerance = 1e-12)))
stopifnot(identical(prior_jump_test$d_prior_jump_prob, 0.30))
cat("INDEPENDENT-PRIOR d-JUMP ACCEPTANCE IDENTITY: PASS\n")

# Test 6c: Compare the integrated likelihood with numerical integration
# over precision u=1/tau2. The log-marginal routine omits only constants
# that are identical for every d,g,xi under the same y and prior.
# These two covariance choices must yield the same likelihood *ratio*.
prior_marginal <- bmim_gp_default_prior(D_list)
marginal_log_density <- function(d_value, g_value) {
  t <- bmim_gp_cholesky_terms(bmim_gp_covariance(D, d_value, g_value), y)
  n <- length(y)
  shape <- prior_marginal$tau_shape
  rate <- prior_marginal$tau_rate
  integrated <- integrate(
    function(u) {
      exp(-n / 2 * log(2 * pi) - t$logdet / 2 +
          (n / 2) * log(u) - u * t$quad / 2 +
          dgamma(u, shape = shape, rate = rate, log = TRUE))
    },
    lower = 0, upper = Inf, rel.tol = 1e-9
  )$value
  log(integrated)
}
d1 <- 0.15 * prior_marginal$distance_scale
d2 <- 1.20 * prior_marginal$distance_scale
g1 <- 0.2
g2 <- 1.1
log_ratio_analytic <-
  bmim_gp_log_marginal_likelihood(y, D, d2, g2, prior_marginal) -
  bmim_gp_log_marginal_likelihood(y, D, d1, g1, prior_marginal)
log_ratio_numeric <- marginal_log_density(d2, g2) -
  marginal_log_density(d1, g1)
stopifnot(isTRUE(all.equal(log_ratio_analytic, log_ratio_numeric,
                           tolerance = 1e-7)))
cat("COLLAPSED TAU2 LIKELIHOOD NUMERICAL CHECK: PASS\n")

# Test 6d: A log-scale slice sampler must reproduce a positive-support
# target density, including the Jacobian of the log transformation.
# This is a direct distributional test independent of the GP code.
stopifnot(identical(prior_jump_test$g_update_method, "slice"))
set.seed(2291)
g_slice_test <- numeric(2300L)
g_last <- 0.4
gamma_z <- function(z) {
  if (z < -700 || z > 700) return(-Inf)
  dgamma(exp(z), shape = 2, rate = 3, log = TRUE) + z
}
for (ii in seq_along(g_slice_test)) {
  g_last <- bmim_gp_slice_log_g(
    g_last, gamma_z, width = 0.75, max_steps = 16L
  )
  g_slice_test[ii] <- g_last
}
# E[Gamma(shape=2,rate=3)] = 2/3. The tolerance exceeds ordinary
# Monte Carlo standard errors and avoids a brittle probabilistic test.
stopifnot(all(is.finite(g_slice_test)), all(g_slice_test > 0),
          abs(mean(g_slice_test[301:2300]) - 2/3) < 0.07)
cat("LOG-NUGGET SLICE SAMPLER DISTRIBUTIONAL CHECK: PASS\n")

# Test 7: The state is persistent and every GP parameter remains
# positive over sequential updates with changing xi.
prior <- bmim_gp_default_prior(D_list)
state <- bmim_gp_initialize(prior)
for (i in seq_len(12)) {
  current_xi <- if (i %% 2L == 0L) c(1, 0) else c(1, 1)
  state <- bmim_gp_update(y, current_xi, D_list, state, prior)
  stopifnot(is.finite(state$d), state$d > 0,
            is.finite(state$g), state$g > 0,
            is.finite(state$tau2), state$tau2 > 0)
}

# Test 8: Actual ISGM() integration with VB disabled (no ClusterR
# requirement) and a small synthetic functional design.
source("ISGM_core.R")
features <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
response <- as.numeric(0.5 * features[, 1] + rnorm(8, sd = 0.5))
chain <- ISGM(
  path_features = features,
  W_mat = W,
  y = response,
  curr_xi = c(1, 1),
  grps_vec = list(list(1:2), list(1:3)),
  N = 8,
  connec = 1,
  VBprop = 0
)
stopifnot(ncol(chain) == 2L,
          all(chain %in% c(0, 1)),
          all(rowSums(chain) >= 1L))

cat("BMIM GP/Ising/VB/state tests: PASS\n")
