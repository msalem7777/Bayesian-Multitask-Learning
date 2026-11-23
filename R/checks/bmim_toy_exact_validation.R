# ============================================================================
# BMIM toy audit: independent enumeration + GP/Ising/VB MCMC + mixing checks
#
# Designed for the corrected R/bmim_gp_updates.R module.
# Use from the R directory:
#   & $rs -e "source('checks/bmim_toy_exact_validation.R')"
#
# This is a STANDALONE CHECK, not a production simulation or a model refit.
# The VB probabilities are held fixed as required by the conditional
# empirical-Bayes target. It does NOT validate the upstream data reuse in VB,
# cross-validation domain selection, or the post-selection score test.
#
# Four independent chains are started in different selection states.
# Independent 2-D numerical integration gives the posterior reference for all
# seven allowed nonempty configurations after integrating out tau2 analytically.
#
# No external R packages are needed.
# ============================================================================

# ----- 0. Load the actual standalone GP update implementation ---------------
script_files <- Filter(
  Negate(is.null),
  lapply(sys.frames(), function(frame) frame$ofile)
)
if (length(script_files) < 1L) {
  stop("Source this test from a real .R file, not via copy-paste expressions.")
}
script_file <- normalizePath(tail(script_files, 1L)[[1L]])
r_dir <- dirname(dirname(script_file))
module_file <- file.path(r_dir, "bmim_gp_updates.R")
if (!file.exists(module_file)) {
  stop("Missing bmim_gp_updates.R at: ", module_file)
}
source(module_file)
rm(script_files, script_file, r_dir, module_file)

# ----- 1. Establish a deterministic small data set --------------------------
set.seed(20261118)
n_subjects <- 9L
n_nodes <- 3L
z <- seq(-1.4, 1.4, length.out = n_subjects)

# Two functional time measurements per node, so X has 9 rows x 6 columns.
X <- cbind(
  z + 0.25 * sin(2 * z),
  0.7 * z + 0.1 * cos(3 * z),
  sin(2.2 * z),
  cos(1.6 * z) + 0.2 * z,
  z^2 + 0.1 * sin(3 * z),
  cos(2.8 * z) - 0.25 * z
)
X_s <- scale(X)
y <- 0.8 * X_s[, 1] - 0.5 * X_s[, 4] +
  0.2 * X_s[, 5] + 0.8 * sin(4.3 * z)
y <- as.numeric(scale(y))

loc_slices <- list(1:2, 3:4, 5:6)
D_by_loc <- bmim_gp_precompute_distances(X_s, loc_slices)
prior <- bmim_gp_default_prior(D_by_loc)

# Example of a network with positive pairwise interactions.
W_mat <- matrix(c(
  0.00, 0.30, 0.15,
  0.30, 0.00, 0.20,
  0.15, 0.20, 0.00
), nrow = n_nodes, byrow = TRUE)

# These are FIXED illustrative VB-derived empirical inclusion probabilities.
# This audit deliberately does not refit VB from the same response.
p_prop <- c(0.67, 0.36, 0.58)

# A model is a binary indicator vector. The actual ISGM algorithm forbids
# the all-zero model; enumerate exactly the 2^3 - 1 allowed configurations.
all_configs <- as.matrix(expand.grid(rep(list(0:1), n_nodes)))
configs <- all_configs[rowSums(all_configs) > 0, , drop = FALSE]
storage.mode(configs) <- "double"
config_names <- apply(configs, 1L, paste0, collapse = "")

stopifnot(nrow(X_s) == n_subjects)
stopifnot(ncol(X_s) == 6L)
stopifnot(length(D_by_loc) == n_nodes)
stopifnot(nrow(configs) == 2^n_nodes - 1L)
stopifnot(all(vapply(D_by_loc, function(D) all(dim(D) == c(9L, 9L)), logical(1))))
cat("TOY DIMENSION CHECK: PASS\n")

# ----- 2. Independent implementation of log-prior and log likelihood --------
# This expression is intentionally independent of
# bmim_gp_selection_log_target(), making it a useful reference.
independent_log_selection <- function(xi, d, g, tau2) {
  s <- 2 * xi - 1
  ising_log <- 0.5 * sum(outer(s, s) * W_mat)
  vb_log <- sum(ifelse(xi == 1, log(p_prop), log1p(-p_prop)))

  D <- matrix(0, nrow = n_subjects, ncol = n_subjects)
  for (p in which(xi == 1)) {
    D <- D + D_by_loc[[p]]
  }
  covariance <- exp(-D / d) + diag(g, n_subjects)
  determinant_log <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
  quadratic_form <- as.numeric(crossprod(y, solve(covariance, y)))

  ising_log + vb_log - 0.5 * determinant_log - 0.5 * quadratic_form / tau2
}

reference_state <- list(d = prior$distance_scale, g = 0.6, tau2 = 1.3)
for (k in seq_len(nrow(configs))) {
  xi <- configs[k, ]
  modular_log <- bmim_gp_selection_log_target(
    xi = xi,
    y = y,
    D_by_loc = D_by_loc,
    state = reference_state,
    W_mat = W_mat,
    p_prop = p_prop
  )
  direct_log <- independent_log_selection(
    xi, reference_state$d, reference_state$g, reference_state$tau2
  )
  stopifnot(isTRUE(all.equal(modular_log, direct_log, tolerance = 1e-9)))
}
cat("INDEPENDENT LOG-TARGET CHECK: PASS\n")

# ----- 3. Detailed balance for EVERY allowable one-coordinate edge ---------
# Two configurations differing in one node have symmetric deterministic
# proposals, so their equilibrium probability flow must match.
for (a in seq_len(nrow(configs))) {
  for (b in seq_len(nrow(configs))) {
    if (sum(abs(configs[a, ] - configs[b, ])) != 1) next
    ua <- independent_log_selection(
      configs[a, ], reference_state$d, reference_state$g, reference_state$tau2
    )
    ub <- independent_log_selection(
      configs[b, ], reference_state$d, reference_state$g, reference_state$tau2
    )
    flow_forward <- ua + min(0, ub - ua)
    flow_reverse <- ub + min(0, ua - ub)
    stopifnot(isTRUE(all.equal(flow_forward, flow_reverse, tolerance = 1e-11)))
  }
}
cat("ALL ALLOWABLE NODE-FLIP DETAILED-BALANCE CHECKS: PASS\n")

# ----- 4. Integrate tau2 analytically, then d and g numerically -------------
# If tau2 ~ IG(a0,b0), then
#   p(y | xi,d,g) proportional to
#     det(C)^(-1/2) * (b0 + y' C^-1 y / 2)^(-(a0 + n/2)).
# This removes Monte Carlo error from the variance integration.
# d and g are integrated on a deterministic midpoint grid in log coordinates.
# Repeating at two resolutions checks whether the grid is sufficiently fine.
log_sum_exp <- function(a) {
  m <- max(a)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(a - m)))
}

normalize_log_weights <- function(logw) {
  exp(logw - log_sum_exp(logw))
}

midpoint_log_grid <- function(low, high, n) {
  width <- (high - low) / n
  low + ((seq_len(n) - 0.5) * width)
}

integrated_reference <- function(grid_size) {
  n_d <- as.integer(grid_size)
  n_g <- as.integer(grid_size)
  log_d_low <- log(prior$distance_scale * 1e-5)
  log_d_high <- log(prior$distance_scale * 12)
  log_g_low <- log(1e-5)
  log_g_high <- log(18)
  dgrid <- midpoint_log_grid(log_d_low, log_d_high, n_d)
  ggrid <- midpoint_log_grid(log_g_low, log_g_high, n_g)
  delta_d <- (log_d_high - log_d_low) / n_d
  delta_g <- (log_g_high - log_g_low) / n_g

  result <- rep(-Inf, nrow(configs))
  a_post <- prior$tau_shape + n_subjects / 2
  distance_by_config <- lapply(seq_len(nrow(configs)), function(k) {
    bmim_gp_selected_distance(D_by_loc, configs[k, ])
  })
  network_vb_log <- vapply(seq_len(nrow(configs)), function(k) {
    xi <- configs[k, ]
    s <- 2 * xi - 1
    0.5 * sum(outer(s, s) * W_mat) +
      sum(ifelse(xi == 1, log(p_prop), log1p(-p_prop)))
  }, numeric(1))

  for (k in seq_len(nrow(configs))) {
    terms <- numeric(n_d * n_g)
    pos <- 0L
    D <- distance_by_config[[k]]
    for (ld in dgrid) {
      d <- exp(ld)
      lp_d <- bmim_gp_log_prior_d(d, prior)
      for (lg in ggrid) {
        g <- exp(lg)
        lp_g <- bmim_gp_log_prior_g(g, prior)
        C <- exp(-D / d) + diag(g, n_subjects)
        R <- chol(C)
        v <- forwardsolve(t(R), y)
        quad <- sum(v^2)
        logdet <- 2 * sum(log(diag(R)))
        pos <- pos + 1L
        # The +ld+lg terms are Jacobians for d=exp(ld), g=exp(lg).
        terms[pos] <- lp_d + lp_g + ld + lg +
          network_vb_log[k] - 0.5 * logdet -
          a_post * log(prior$tau_rate + quad / 2)
      }
    }
    result[k] <- log_sum_exp(terms) + log(delta_d) + log(delta_g)
  }
  normalize_log_weights(result)
}

cat("Computing independent low/high-resolution numerical reference...\n")
reference_low <- integrated_reference(grid_size = 40L)
reference_high <- integrated_reference(grid_size = 64L)
quadrature_difference <- max(abs(reference_high - reference_low))
cat(sprintf("Largest grid-refinement probability change: %.5f\n", quadrature_difference))
stopifnot(quadrature_difference < 0.025)
exact_reference <- reference_high
names(exact_reference) <- config_names
cat("INDEPENDENT NUMERICAL REFERENCE: PASS\n")

# ----- 5. Run actual modular GP + deterministic node-flip transitions -------
# Each chain starts in a different configuration. The GP state persists
# across sweeps. The full node target includes both Ising and VB factors.
run_chain <- function(initial_xi, chain_number, n_iteration, burnin) {
  xi <- as.numeric(initial_xi)
  state <- bmim_gp_initialize(prior)
  initial_d_multiplier <- c(0.1, 0.5, 2.0, 5.0)[chain_number]
  initial_g <- c(0.07, 0.25, 1.5, 4.0)[chain_number]
  initial_tau2 <- c(0.25, 0.75, 2.0, 6.0)[chain_number]
  state$d <- prior$distance_scale * initial_d_multiplier
  state$g <- initial_g
  state$tau2 <- initial_tau2

  kept <- n_iteration - burnin
  xi_draws <- matrix(NA_real_, nrow = kept, ncol = n_nodes)
  log_parameter_draws <- matrix(NA_real_, nrow = kept, ncol = 3L)
  colnames(log_parameter_draws) <- c("log_d", "log_g", "log_tau2")
  accepted_node <- 0L
  attempted_node <- 0L
  accepted_d <- 0L
  accepted_g <- 0L

  for (iteration in seq_len(n_iteration)) {
    state <- bmim_gp_update(
      y = y, xi = xi, D_by_loc = D_by_loc, state = state, prior = prior
    )
    accepted_d <- accepted_d + as.integer(state$accepted_d)
    accepted_g <- accepted_g + as.integer(state$accepted_g)

    # Exactly the same deterministic coordinate flip and acceptance formula
    # as the ISGM() loop, with all GP parameters held fixed.
    for (j in seq_len(n_nodes)) {
      proposed <- xi
      proposed[j] <- 1 - proposed[j]
      if (sum(proposed) == 0) next
      attempted_node <- attempted_node + 1L
      old_log <- bmim_gp_selection_log_target(
        xi, y, D_by_loc, state, W_mat, p_prop
      )
      new_log <- bmim_gp_selection_log_target(
        proposed, y, D_by_loc, state, W_mat, p_prop
      )
      if (log(runif(1)) < min(0, new_log - old_log)) {
        xi <- proposed
        accepted_node <- accepted_node + 1L
      }
    }

    if (iteration > burnin) {
      idx <- iteration - burnin
      xi_draws[idx, ] <- xi
      log_parameter_draws[idx, ] <- log(c(state$d, state$g, state$tau2))
    }
  }
  list(
    xi = xi_draws,
    gp = log_parameter_draws,
    node_acceptance = accepted_node / max(1, attempted_node),
    d_acceptance = accepted_d / n_iteration,
    g_acceptance = accepted_g / n_iteration # slice draw valid, not MH acceptance
  )
}

n_chains <- 4L
n_iteration <- 6500L
burnin <- 1500L
initial_states <- list(c(1,0,0), c(0,1,0), c(0,0,1), c(1,1,1))
cat(sprintf("Running %d chains, %d sweeps each (%d burn-in)...\n",
            n_chains, n_iteration, burnin))
chains <- lapply(seq_len(n_chains), function(i) {
  run_chain(initial_states[[i]], i, n_iteration, burnin)
})

# ----- 6. Check exact model probabilities against empirical frequencies -----
all_xi <- do.call(rbind, lapply(chains, function(x) x$xi))
visited_configs <- apply(all_xi, 1L, paste0, collapse = "")
counts <- table(factor(visited_configs, levels = config_names))
observed <- as.numeric(counts) / nrow(all_xi)
names(observed) <- config_names
results <- data.frame(
  configuration = config_names,
  oracle_probability = round(as.numeric(exact_reference), 5),
  mcmc_probability = round(as.numeric(observed), 5),
  absolute_error = round(abs(observed - exact_reference), 5),
  row.names = NULL
)
cat("\nALL NONEMPTY NODE CONFIGURATIONS:\n")
print(results, row.names = FALSE)
max_configuration_error <- max(abs(observed - exact_reference))
cat(sprintf("Maximum probability discrepancy: %.5f\n", max_configuration_error))

# ----- 7. Four-chain split-Rhat and autocorrelation ESS diagnostics ----------
# Split Rhat compares chain means with within-chain variability.
# ESS = effective sample size: an approximation of the number of independent
# draws represented by correlated draws. These checks are diagnostics, not
# mathematical proofs of convergence.
compute_split_rhat <- function(samples_by_chain) {
  half <- floor(length(samples_by_chain[[1L]]) / 2)
  split_chains <- unlist(lapply(samples_by_chain, function(x) {
    list(x[seq_len(half)], x[half + seq_len(half)])
  }), recursive = FALSE)
  chain_means <- vapply(split_chains, mean, numeric(1))
  chain_vars <- vapply(split_chains, var, numeric(1))
  within <- mean(chain_vars)
  between <- half * var(chain_means)
  if (within == 0) {
    if (between == 0) return(NA_real_)  # no movement; not validated
    return(Inf)
  }
  variance_estimate <- (half - 1) / half * within + between / half
  sqrt(variance_estimate / within)
}

compute_ess <- function(samples_by_chain, max_lag = 240L) {
  # Autocorrelation-based Geyer initial-positive-pair estimate.
  if (any(vapply(samples_by_chain, function(x) var(x) < 1e-14, logical(1)))) {
    return(0)
  }
  n <- length(samples_by_chain[[1L]])
  lagmax <- min(as.integer(max_lag), n - 1L)
  correlations <- vapply(samples_by_chain, function(x) {
    as.numeric(stats::acf(x, lag.max = lagmax, plot = FALSE)$acf)[-1L]
  }, numeric(lagmax))
  if (is.null(dim(correlations))) {
    correlations <- matrix(correlations, ncol = length(samples_by_chain))
  }
  mean_rho <- rowMeans(correlations)
  n_pairs <- floor(length(mean_rho) / 2)
  sum_positive <- 0
  for (k in seq_len(n_pairs)) {
    pair_sum <- mean_rho[2L * k - 1L] + mean_rho[2L * k]
    if (!is.finite(pair_sum) || pair_sum <= 0) break
    sum_positive <- sum_positive + pair_sum
  }
  tau_int <- max(1, 1 + 2 * sum_positive)
  min(n * length(samples_by_chain), n * length(samples_by_chain) / tau_int)
}

parameters <- c("xi_1", "xi_2", "xi_3", "log_d", "log_g", "log_tau2")
monitor <- lapply(seq_along(parameters), function(k) {
  sample_list <- lapply(chains, function(ch) {
    if (k <= n_nodes) ch$xi[, k] else ch$gp[, k - n_nodes]
  })
  c(rhat = compute_split_rhat(sample_list), ess = compute_ess(sample_list))
})
monitor <- do.call(rbind, monitor)
diagnostics <- data.frame(
  quantity = parameters,
  split_rhat = round(monitor[, "rhat"], 4),
  effective_sample_size = round(monitor[, "ess"], 1),
  row.names = NULL
)
cat("\nCONVERGENCE DIAGNOSTICS (base-R estimates):\n")
print(diagnostics, row.names = FALSE)
cat("\nMCMC UPDATE DIAGNOSTICS (not convergence measures):\n")
print(data.frame(
  chain = seq_len(n_chains),
  node_flip_acceptance = round(vapply(chains, `[[`, numeric(1), "node_acceptance"), 3),
  d_acceptance = round(vapply(chains, `[[`, numeric(1), "d_acceptance"), 3),
  g_slice_update = round(vapply(chains, `[[`, numeric(1), "g_acceptance"), 3)
), row.names = FALSE)

# ----- 8. Make failure informative rather than reporting a false PASS --------
# These are illustrative smoke-test thresholds, NOT guarantees of calibration.
conditions <- c(
  quadrature_stable = quadrature_difference < 0.025,
  target_frequencies_match = max_configuration_error < 0.06,
  all_split_rhat_below_1_05 = all(is.finite(monitor[, "rhat"])) &&
    all(monitor[, "rhat"] < 1.05),
  effective_sample_size_above_250 = all(is.finite(monitor[, "ess"])) &&
    all(monitor[, "ess"] > 250)
)
cat("\nVALIDATION CHECKS:\n")
print(conditions)
if (!all(conditions)) {
  stop("TOY CHECKS NEED REVIEW: see frequencies and mixing diagnostics above.")
}
cat("\nBMIM TOY: TARGET + FOUR-CHAIN EXPLORATION CHECKS PASS\n")
