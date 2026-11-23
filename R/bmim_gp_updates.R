# ============================================================
# bmim_gp_updates.R
# GP-parameter updates with analytically marginalized scale for d,g.
#
# Model targeted by this module (conditional on preprocessed y):
#   y | xi, d, g, tau2 ~ N(0, tau2 * [K(xi,d) + g*I])
#   K_ij(xi,d) = exp(-sum_p xi_p * D_p[i,j] / d)
#   tau2 ~ Inverse-Gamma(alpha_tau/2, q_tau/2)
#   d / distance_scale ~ 0.5*Gamma(1,20) + 0.5*Gamma(10,10)
#   g ~ Gamma(1,1)
#
# IMPORTANT: This is a conditional empirical-Bayes model when y
# and the VB probabilities have been calculated from the data.
# It is NOT a sampler for the original paper's entire beta/beta0/W
# hierarchy. All functions use only base R.
# ============================================================

# Precompute, for each functional location, the matrix of squared
# distances between subjects. Every retained time column receives
# a factor 1 / number_of_retained_time_columns_for_location.
bmim_gp_precompute_distances <- function(X_s, loc_slices) {
  X_s <- as.matrix(X_s)
  storage.mode(X_s) <- "double"
  if (!is.numeric(X_s) || any(!is.finite(X_s))) {
    stop("X_s must contain only finite standardized numeric values.")
  }
  if (nrow(X_s) < 2L || length(loc_slices) < 1L) {
    stop("Need at least two subjects and one functional location.")
  }
  n <- nrow(X_s)
  out <- vector("list", length(loc_slices))
  for (p in seq_along(loc_slices)) {
    ids <- loc_slices[[p]]
    if (length(ids) == 0L) {
      out[[p]] <- matrix(0, n, n)
      next
    }
    # Scale each node by its number of retained time columns.
    Xp <- X_s[, ids, drop = FALSE] / sqrt(length(ids))

    # Squared Euclidean distance between subjects.
    squared_norm <- rowSums(Xp * Xp)
    D <- outer(squared_norm, squared_norm, "+") -
      2 * tcrossprod(Xp)

    # Enforce exact symmetry while preserving matrix dimensions.
    D <- (D + t(D)) / 2

    # Remove tiny negative values caused by floating-point rounding.
    # Matrix indexing preserves the dimensions of D.
    D[D < 0] <- 0

    # Each subject has zero distance from itself.
    diag(D) <- 0

    out[[p]] <- D
  }
  out
}

# Explicit prior choices. The distance scale is fixed using inputs only,
# not selected using the response or changed during the chain.
# For a positive squared-distance matrix, it is the median positive
# distance under the all-nodes configuration.
bmim_gp_default_prior <- function(D_by_loc,
                                  alpha_tau = 4,
                                  q_tau = 2,
                                  initial_g = 0.1,
                                  proposal_sd_d = 0.25,
                                  proposal_sd_g = 0.20,
                                  d_prior_jump_prob = 0.30,
                                  update_g = TRUE,
                                  g_update_method = "slice",
                                  g_slice_width_log = 0.75,
                                  g_slice_max_steps = 16L,
                                  g_slice_max_shrink = 200L,
                                  null_sigma_shape = alpha_tau / 2,
                                  null_sigma_rate = q_tau / 2) {
  D_full <- Reduce("+", D_by_loc)
  values <- D_full[upper.tri(D_full)]
  values <- values[is.finite(values) & values > 0]
  distance_scale <- if (length(values) > 0L) median(values) else 1
  if (!is.finite(distance_scale) || distance_scale <= 0) {
    stop("Invalid distance scale.")
  }
  if (alpha_tau <= 0 || q_tau <= 0 || initial_g <= 0 ||
      proposal_sd_d <= 0 || proposal_sd_g <= 0 ||
      !is.finite(d_prior_jump_prob) ||
      d_prior_jump_prob < 0 || d_prior_jump_prob > 1 ||
      length(g_update_method) != 1L ||
      !(g_update_method %in% c("slice", "rw")) ||
      !is.finite(g_slice_width_log) || g_slice_width_log <= 0 ||
      length(g_slice_max_steps) != 1L || !is.finite(g_slice_max_steps) ||
      g_slice_max_steps < 1L || g_slice_max_steps != as.integer(g_slice_max_steps) ||
      length(g_slice_max_shrink) != 1L || !is.finite(g_slice_max_shrink) ||
      g_slice_max_shrink < 1L || g_slice_max_shrink != as.integer(g_slice_max_shrink) ||
      !is.finite(null_sigma_shape) || null_sigma_shape <= 0 ||
      !is.finite(null_sigma_rate) || null_sigma_rate <= 0) {
    stop("GP prior and proposal constants must be positive.")
  }
  list(
    tau_shape = alpha_tau / 2,
    tau_rate = q_tau / 2,
    distance_scale = distance_scale,
    d_mix_weights = c(0.5, 0.5),
    d_mix_shapes = c(1, 10),
    d_mix_rates = c(20, 10),
    g_shape = 1,
    g_rate = 1,
    initial_g = initial_g,
    proposal_sd_d = proposal_sd_d,
    proposal_sd_g = proposal_sd_g,
    d_prior_jump_prob = d_prior_jump_prob,
    update_g = update_g,
    g_update_method = g_update_method,
    g_slice_width_log = g_slice_width_log,
    g_slice_max_steps = as.integer(g_slice_max_steps),
    g_slice_max_shrink = as.integer(g_slice_max_shrink),
    null_sigma_shape = null_sigma_shape,
    null_sigma_rate = null_sigma_rate
  )
}

# Sum the per-location squared distances for selected locations.
bmim_gp_selected_distance <- function(D_by_loc, xi) {
  if (length(D_by_loc) != length(xi) ||
      any(!xi %in% c(0, 1))) {
    stop("xi must be binary, with one entry per location.")
  }
  n <- nrow(D_by_loc[[1L]])
  D <- matrix(0, n, n)
  for (p in which(xi == 1)) D <- D + D_by_loc[[p]]
  D
}

# The covariance used everywhere, including GP updates and node flips.
# For the empty selection, exp(0) gives a constant kernel, but ISGM()
# may separately forbid empty selections as part of its model support.
bmim_gp_covariance <- function(D, d, g) {
  if (length(d) != 1L || length(g) != 1L ||
      !is.finite(d) || !is.finite(g) || d <= 0 || g <= 0) {
    stop("d and g must be positive finite scalars.")
  }
  K <- exp(-D / d)
  diag(K) <- diag(K) + g
  K
}

# Gaussian quadratic form and log determinant using one Cholesky factor.
bmim_gp_cholesky_terms <- function(K, y) {
  y <- as.numeric(y)
  if (!is.matrix(K) || nrow(K) != ncol(K) ||
      nrow(K) != length(y)) stop("K and y dimensions disagree.")
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  z <- forwardsolve(t(R), y)
  quad <- sum(z * z)
  logdet <- 2 * sum(log(diag(R)))
  if (!is.finite(quad) || !is.finite(logdet)) return(NULL)
  list(quad = quad, logdet = logdet)
}

# Log Gaussian likelihood up to constants independent of the node
# selection and GP hyperparameters for fixed tau2.
bmim_gp_log_likelihood <- function(y, D, d, g, tau2) {
  if (!is.finite(tau2) || tau2 <= 0) return(-Inf)
  K <- bmim_gp_covariance(D, d, g)
  terms <- bmim_gp_cholesky_terms(K, y)
  if (is.null(terms)) return(-Inf)
  -0.5 * (terms$logdet + terms$quad / tau2)
}

# Marginal Gaussian likelihood after analytically integrating tau2.
#
# If y | C,tau2 ~ N(0,tau2*C) and tau2 ~ IG(a,b) with
# density proportional to tau2^(-a-1)*exp(-b/tau2), then
# p(y | C) is proportional to
#   |C|^(-1/2) * (b + 0.5*y^T C^(-1) y)^(-(a+n/2)).
# Constants independent of d,g,xi are omitted.
# IMPORTANT: The d,g MH updates use this marginalized likelihood.
# The node flips still use the existing conditional likelihood given
# tau2; the order is: collapsed d,g -> Gibbs tau2 -> node flips.
bmim_gp_log_marginal_likelihood <- function(y, D, d, g, prior) {
  K <- bmim_gp_covariance(D, d, g)
  terms <- bmim_gp_cholesky_terms(K, y)
  if (is.null(terms)) return(-Inf)
  a <- prior$tau_shape
  b <- prior$tau_rate
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("Inverse-gamma prior shape and rate must be positive.")
  }
  -0.5 * terms$logdet -
    (a + length(y) / 2) * log(b + terms$quad / 2)
}

# The tgp-style mixture Gamma is imposed on the dimensionless quantity
# d / distance_scale; this scaling is an EXPLICIT new modeling choice.
bmim_gp_log_prior_d <- function(d, prior) {
  if (!is.finite(d) || d <= 0) return(-Inf)
  x <- d / prior$distance_scale
  log_terms <- log(prior$d_mix_weights) +
    dgamma(x, shape = prior$d_mix_shapes,
           rate = prior$d_mix_rates, log = TRUE) -
    log(prior$distance_scale)
  m <- max(log_terms)
  m + log(sum(exp(log_terms - m)))
}

bmim_gp_log_prior_g <- function(g, prior) {
  if (!is.finite(g) || g <= 0) return(-Inf)
  dgamma(g, shape = prior$g_shape, rate = prior$g_rate,
         log = TRUE)
}

# Metropolis proposal on log(parameter). The log(proposed/current)
# term is the Jacobian correction from changing parameterization.
bmim_gp_log_mh_ratio <- function(logpost_current, logpost_proposed,
                                 current, proposed) {
  if (!is.finite(current) || !is.finite(proposed) ||
      current <= 0 || proposed <= 0) stop("Positive states required.")
  logpost_proposed - logpost_current + log(proposed) - log(current)
}

bmim_gp_initialize <- function(prior) {
  list(d = prior$distance_scale,
       g = prior$initial_g,
       tau2 = 1,
       accepted_d = FALSE,
       accepted_g = FALSE)
}

# Draw g from its one-dimensional conditional with a slice update.
# A slice is the region where the log-posterior exceeds a randomly chosen
# horizontal level. Stepping out and then shrinking that region avoids
# relying on a tiny log-scale random walk. The calculation is performed
# on z=log(g), so log p_z(z)=log p_g(exp(z)) + z (Jacobian).
# This is a standard fixed-width, bounded-stepping-out slice kernel;
# the random allocation of left/right expansions preserves reversibility.
bmim_gp_slice_log_g <- function(g, log_target_z, width = 0.75,
                                max_steps = 16L, max_shrink = 200L) {
  if (length(g) != 1L || !is.finite(g) || g <= 0 ||
      length(width) != 1L || !is.finite(width) || width <= 0 ||
      length(max_steps) != 1L || max_steps < 1L ||
      length(max_shrink) != 1L || max_shrink < 1L) {
    stop("Invalid log-g slice-sampler settings.")
  }
  z0 <- log(g)
  l0 <- log_target_z(z0)
  if (!is.finite(l0)) stop("Current nugget has nonfinite log posterior.")
  level <- l0 - rexp(1L)
  left <- z0 - runif(1L) * width
  right <- left + width
  j <- sample.int(as.integer(max_steps), 1L) - 1L
  k <- as.integer(max_steps) - 1L - j
  while (j > 0L && log_target_z(left) > level) {
    left <- left - width
    j <- j - 1L
  }
  while (k > 0L && log_target_z(right) > level) {
    right <- right + width
    k <- k - 1L
  }
  for (iteration in seq_len(as.integer(max_shrink))) {
    candidate <- runif(1L, min = left, max = right)
    candidate_lp <- log_target_z(candidate)
    if (is.finite(candidate_lp) && candidate_lp >= level) {
      return(exp(candidate))
    }
    if (candidate < z0) left <- candidate else right <- candidate
  }
  stop("Nugget slice sampler exhausted max_shrink; increase its limit.")
}

# Posterior log-density on z=log(g), using the scale-integrated likelihood.
# The +z term is essential: it converts a density in g to one in log(g).
bmim_gp_log_g_target_z <- function(z, y, D, d, prior) {
  if (length(z) != 1L || !is.finite(z) || z > 700 || z < -700) {
    return(-Inf)
  }
  g <- exp(z)
  bmim_gp_log_marginal_likelihood(y, D, d, g, prior) +
    bmim_gp_log_prior_g(g, prior) + z
}

# A partially collapsed GP sweep with persistent parameter state.
# 1) Update d from its posterior with tau2 integrated out.
# 2) Update g with a log-scale slice draw from its collapsed posterior.
# 3) Draw tau2 from its exact IG conditional given the updated d,g.
# Then, in ISGM(), update xi conditional on this freshly sampled tau2.
# This ordering preserves the original joint target. Do not insert a
# tau2-conditional d/g move between steps 1 and 3.
bmim_gp_update <- function(y, xi, D_by_loc, state, prior) {
  D <- bmim_gp_selected_distance(D_by_loc, xi)
  d <- state$d
  g <- state$g
  tau2 <- state$tau2
  if (!is.finite(tau2) || tau2 <= 0) stop("Bad prior GP state.")

  # A mixture of two individually valid MH transition kernels:
  # (i) the original small log-random-walk move and
  # (ii) an independent proposal drawn from the two-component prior.
  # The mixing probability is fixed, not adapted to the current state.
  use_prior_jump <- runif(1L) < prior$d_prior_jump_prob
  if (use_prior_jump) {
    component <- sample.int(
      length(prior$d_mix_weights), size = 1L,
      prob = prior$d_mix_weights
    )
    proposed_d <- prior$distance_scale * rgamma(
      1L, shape = prior$d_mix_shapes[component],
      rate = prior$d_mix_rates[component]
    )
    # For an independence proposal q(d*) = p(d*), the proposal-density
    # ratio cancels the identical d-prior ratio in the MH acceptance.
    current <- bmim_gp_log_marginal_likelihood(y, D, d, g, prior)
    proposed <- bmim_gp_log_marginal_likelihood(
      y, D, proposed_d, g, prior
    )
    log_ratio <- proposed - current
  } else {
    current <- bmim_gp_log_marginal_likelihood(y, D, d, g, prior) +
      bmim_gp_log_prior_d(d, prior)
    proposed_d <- exp(log(d) + rnorm(1L, sd = prior$proposal_sd_d))
    proposed <- bmim_gp_log_marginal_likelihood(
      y, D, proposed_d, g, prior
    ) + bmim_gp_log_prior_d(proposed_d, prior)
    log_ratio <- bmim_gp_log_mh_ratio(
      current, proposed, d, proposed_d
    )
  }
  accept_d <- FALSE
  if (is.finite(proposed) && is.finite(current) &&
      is.finite(log_ratio)) {
    accept_d <- log(runif(1L)) < min(0, log_ratio)
    if (accept_d) d <- proposed_d
  }

  # Update g conditional on the accepted d, with tau2 integrated out.
  # Default: univariate slice sampling on log(g), avoiding the tiny accepted
  # random-walk moves found in the initial four-chain diagnostic.
  # The prior's "rw" option is kept as a controlled comparison.
  accept_g <- FALSE
  if (isTRUE(prior$update_g)) {
    method <- if (is.null(prior$g_update_method)) "rw" else prior$g_update_method
    if (identical(method, "slice")) {
      previous_g <- g
      target_z <- function(z) bmim_gp_log_g_target_z(z, y, D, d, prior)
      g <- bmim_gp_slice_log_g(
        g, target_z,
        width = prior$g_slice_width_log,
        max_steps = prior$g_slice_max_steps,
        max_shrink = prior$g_slice_max_shrink
      )
      # Here accepted_g indicates that a valid slice draw was produced,
      # NOT a Metropolis-Hastings acceptance probability.
      accept_g <- is.finite(g) && g > 0 && g != previous_g
    } else if (identical(method, "rw")) {
      current <- bmim_gp_log_marginal_likelihood(y, D, d, g, prior) +
        bmim_gp_log_prior_g(g, prior)
      proposed_g <- exp(log(g) + rnorm(1L, sd = prior$proposal_sd_g))
      proposed <- bmim_gp_log_marginal_likelihood(y, D, d, proposed_g, prior) +
        bmim_gp_log_prior_g(proposed_g, prior)
      if (is.finite(proposed) && is.finite(current)) {
        log_ratio <- bmim_gp_log_mh_ratio(current, proposed, g, proposed_g)
        accept_g <- log(runif(1L)) < min(0, log_ratio)
        if (accept_g) g <- proposed_g
      }
    } else {
      stop("Unknown g update method; use 'slice' or 'rw'.")
    }
  }

  # Draw the GP scale after the collapsed d/g updates.
  terms <- bmim_gp_cholesky_terms(bmim_gp_covariance(D, d, g), y)
  if (is.null(terms)) stop("Current GP covariance is not positive definite.")
  shape_post <- prior$tau_shape + length(y) / 2
  rate_post <- prior$tau_rate + terms$quad / 2
  tau2 <- 1 / rgamma(1L, shape = shape_post, rate = rate_post)

  list(d = d, g = g, tau2 = tau2,
       accepted_d = accept_d, accepted_g = accept_g)
}

# Conditional GP-Ising posterior WITH the explicitly requested
# empirical-Bayes VB factor. The VB probabilities must be computed
# ONCE from the observed data and then held fixed for this chain.
bmim_gp_selection_log_target <- function(xi, y, D_by_loc,
                                         state, W_mat, p_prop) {
  if (length(xi) != length(p_prop) ||
      any(!xi %in% c(0, 1)) ||
      !all(dim(W_mat) == c(length(xi), length(xi)))) {
    stop("Node indicator, VB probability, and network dimensions disagree.")
  }
  if (any(!is.finite(p_prop)) ||
      any(p_prop <= 0 | p_prop >= 1)) {
    stop("All VB probabilities must lie strictly between 0 and 1.")
  }
  spins <- 2 * xi - 1
  ising <- 0.5 * as.numeric(crossprod(spins, W_mat %*% spins))
  log_vb_prior <- sum(dbinom(xi, size = 1, prob = p_prop,
                             log = TRUE))
  D <- bmim_gp_selected_distance(D_by_loc, xi)
  ising + log_vb_prior +
    bmim_gp_log_likelihood(y, D, state$d, state$g, state$tau2)
}
