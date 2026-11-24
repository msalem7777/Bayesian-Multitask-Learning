# BMIM collapsed joint node/domain/global-Bernoulli MCMC driver.
# Requires the two domain modules and validated GP / global modules.
# Source bmim_joint_domains.R to load the complete workflow.

bmim_joint_require_modules <- function(global_indicator = FALSE) {
  required <- c("bmim_gp_precompute_distances", "bmim_gp_default_prior",
                "bmim_gp_initialize", "bmim_gp_update",
                "bmim_gp_log_marginal_likelihood",
                "bmim_gp_covariance", "bmim_gp_cholesky_terms",
                "bmim_scalar_setup", "bmim_scalar_alt_terms",
                "bmim_scalar_log_evidence", "bmim_scalar_draw_beta")
  if (global_indicator) {
    required <- c(required, "bmim_global_log_integrated_tau",
                  "bmim_global_null_logevidence")
  }
  missing <- required[!vapply(required, exists, logical(1L),
                               mode = "function", inherits = TRUE)]
  if (length(missing)) {
    stop("Source prerequisite modules first. Missing: ",
         paste(missing, collapse = ", "))
  }
}

# Define the domain prior for ONE node with T observations.
# Every allowed length has equal prior probability, by default.
# Conditional on length l, every start in {1,...,T-l+1} is equally likely.

bmim_joint_draw_ig <- function(shape, rate) {
  1 / rgamma(1L, shape = shape, rate = rate)
}

bmim_joint_sample_d_prior <- function(prior) {
  j <- sample.int(length(prior$d_mix_weights), 1L,
                  prob = prior$d_mix_weights)
  prior$distance_scale * rgamma(1L, shape = prior$d_mix_shapes[j],
                                rate = prior$d_mix_rates[j])
}

# Draw GP and null scales AFTER every collapsed xi/A/d/g/b update.
# This restores the full joint state without conditioning collapsed moves
# on a stale tau2. The non-active model's scale is an auxiliary prior draw.
bmim_joint_draw_scales <- function(b, y, D_total, gp, prior,
                                    scalar = NULL, beta = NULL) {
  if (b == 1L) {
    terms <- bmim_scalar_alt_terms(y, D_total, gp$d, gp$g, scalar, beta)
    if (is.null(terms)) stop("Invalid covariance in joint-domain scale draw.")
    tau2 <- bmim_joint_draw_ig(
      prior$tau_shape + terms$dimension / 2,
      prior$tau_rate + terms$quad / 2
    )
    sigma2 <- bmim_joint_draw_ig(prior$null_sigma_shape,
                                 prior$null_sigma_rate)
  } else {
    terms <- bmim_scalar_null_terms(y, scalar, beta)
    if (is.null(terms)) stop("Invalid covariance in null scale draw.")
    tau2 <- bmim_joint_draw_ig(prior$tau_shape, prior$tau_rate)
    sigma2 <- bmim_joint_draw_ig(
      prior$null_sigma_shape + terms$dimension / 2,
      prior$null_sigma_rate + terms$quad / 2
    )
  }
  list(tau2 = tau2, sigma2 = sigma2)
}

# Main self-contained MCMC on already-standardized full trajectories.
# y_tilde is the historical argument name: with Z supplied, pass the
# original (or fixed-scale) response, NOT regression residuals.
# Z is the exact numeric scalar design, including any desired intercept.
# pi_b=NULL : the GP is always present (joint node/domain sampler).
# pi_b in (0,1) : include the noise-only null, NO node shortlist needed.
# fixed_gp=list(d=...,g=...) : hold these two parameters fixed for an
# independent, discrete-state posterior toy validation; NOT production.
bmim_joint_domain_mcmc <- function(X_s, y_tilde, loc_slices,
                                   W_mat, p_prop, max_domain_length,
                                   n_iter, burnin, n_chains = 4L,
                                   min_domain_length = 1L,
                                   pi_b = NULL,
                                   n_domain_updates_per_sweep = 1L,
                                   domain_prior_jump_prob = 0.50,
                                   gp_prior = NULL,
                                   fixed_gp = NULL,
                                   seed = 20261122L,
                                   Z = NULL,
                                   beta_mode = c("integrated", "sampled"),
                                   beta_prior_mean = NULL,
                                   beta_prior_cov = NULL,
                                   save_beta = FALSE) {
  global_indicator <- !is.null(pi_b)
  bmim_joint_require_modules(global_indicator)
  X_s <- as.matrix(X_s)
  storage.mode(X_s) <- "double"
  y <- as.numeric(y_tilde)
  beta_mode <- match.arg(beta_mode)
  if (!is.logical(save_beta) || length(save_beta) != 1L || is.na(save_beta))
    stop("save_beta must be TRUE or FALSE.")
  P <- length(loc_slices)
  n <- nrow(X_s)
  is_pos_int <- function(z) length(z) == 1L && is.numeric(z) &&
    is.finite(z) && z >= 1L && z == as.integer(z)
  if (n < 2L || P < 1L || ncol(X_s) < P ||
      any(!is.finite(X_s)) || length(y) != n || any(!is.finite(y)) ||
      any(vapply(loc_slices, length, integer(1L)) == 0L) ||
      any(unlist(loc_slices) < 1L) ||
      any(unlist(loc_slices) > ncol(X_s)) ||
      !identical(sort(as.integer(unlist(loc_slices))), seq_len(ncol(X_s))) ||
      !all(dim(W_mat) == c(P,P)) || any(!is.finite(W_mat)) ||
      max(abs(W_mat - t(W_mat))) > 1e-9 ||
      length(p_prop) != P || any(!is.finite(p_prop)) ||
      any(p_prop <= 0 | p_prop >= 1) ||
      !is_pos_int(max_domain_length) || !is_pos_int(min_domain_length) ||
      !is_pos_int(n_iter) || n_iter < 2L ||
      length(burnin) != 1L || !is.numeric(burnin) ||
      !is.finite(burnin) || burnin < 0 || burnin >= n_iter ||
      burnin != as.integer(burnin) || !is_pos_int(n_chains) ||
      !is_pos_int(n_domain_updates_per_sweep) ||
      length(domain_prior_jump_prob) != 1L ||
      !is.finite(domain_prior_jump_prob) ||
      domain_prior_jump_prob < 0 || domain_prior_jump_prob > 1 ||
      !is_pos_int(seed)) {
    stop("Invalid joint-domain dimensions, priors, or MCMC controls.")
  }
  if (global_indicator && (length(pi_b) != 1L ||
      !is.finite(pi_b) || pi_b <= 0 || pi_b >= 1)) {
    stop("pi_b must be in (0,1), or NULL to fix b=1.")
  }
  if (!is.null(fixed_gp) && (is.null(fixed_gp$d) ||
      is.null(fixed_gp$g) || length(fixed_gp$d) != 1L ||
      length(fixed_gp$g) != 1L || !is.finite(fixed_gp$d) ||
      !is.finite(fixed_gp$g) || fixed_gp$d <= 0 || fixed_gp$g <= 0)) {
    stop("fixed_gp requires positive finite d and g.")
  }
  scalar <- bmim_scalar_setup(Z, n, beta_prior_mean, beta_prior_cov)
  sampling_beta <- !is.null(scalar) && beta_mode == "sampled"
  recording_beta <- !is.null(scalar) && (sampling_beta || save_beta)
  domain_priors <- lapply(loc_slices, function(indices) {
    bmim_domain_prior(length(indices), max_domain_length,
                      min_domain_length)
  })
  # The length-scale prior's distance scale is based ONLY on complete inputs,
  # is computed ONCE, and does not drift as intervals move.
  if (is.null(gp_prior)) {
    complete_distances <- bmim_gp_precompute_distances(X_s, loc_slices)
    gp_prior <- bmim_gp_default_prior(complete_distances)
  }
  if (global_indicator && (is.null(gp_prior$null_sigma_shape) ||
      is.null(gp_prior$null_sigma_rate))) {
    stop("Global model requires a proper independent null variance prior.")
  }
  # Integrated null evidence is constant across all MCMC states.
  # In sampled mode the shared beta changes, so recompute each sweep.
  log_null <- if (global_indicator && !sampling_beta) {
    if (is.null(scalar)) bmim_global_null_logevidence(y, gp_prior) else
      bmim_scalar_log_evidence(bmim_scalar_null_terms(y, scalar),
        gp_prior$null_sigma_shape, gp_prior$null_sigma_rate)
  } else NA_real_

  # Local RNG seed; avoid altering the caller's random state.
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (has_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  chains <- vector("list", as.integer(n_chains))
  for (ch in seq_len(as.integer(n_chains))) {
    # Initialization need not be from the exact Ising distribution;
    # MCMC proposals nevertheless use its exact target.
    xi <- as.integer(rbinom(P, 1L, p_prop))
    if (sum(xi) == 0L) xi[which.max(p_prop)] <- 1L
    A <- matrix(NA_integer_, P, 2L, dimnames = list(NULL, c("start", "end")))
    D_by_loc <- vector("list", P)
    for (p in seq_len(P)) {
      A[p, ] <- bmim_domain_draw(domain_priors[[p]])
      D_by_loc[[p]] <- bmim_domain_distance(X_s, loc_slices[[p]], A[p, ])
    }
    D_total <- bmim_gp_selected_distance(D_by_loc, xi)
    gp <- bmim_gp_initialize(gp_prior)
    if (!is.null(fixed_gp)) {
      gp$d <- fixed_gp$d
      gp$g <- fixed_gp$g
    }
    b <- if (global_indicator) as.integer(ch %% 2L == 0L) else 1L
    # beta is held fixed throughout the variance-collapsed sweep in sampled
    # mode. In integrated mode NULL signals analytical marginalization.
    beta <- if (sampling_beta) scalar$mu else NULL
    beta_samples <- if (recording_beta) matrix(NA_real_, n_iter, scalar$s,
      dimnames=list(NULL, scalar$names)) else NULL
    log_marginal <- if (is.null(scalar)) bmim_gp_log_marginal_likelihood else
      function(y, D, d, g, prior) {
        bmim_scalar_log_evidence(bmim_scalar_alt_terms(y,D,d,g,scalar,beta),
                                 prior$tau_shape, prior$tau_rate)
      }
    log_alt_kernel <- function(y, D_total, gp, gp_prior) {
      log_marginal(y, D_total, gp$d, gp$g, gp_prior)
    }
    xi_samples <- matrix(NA_integer_, n_iter, P)
    starts <- matrix(NA_integer_, n_iter, P)
    ends <- matrix(NA_integer_, n_iter, P)
    parameters <- matrix(NA_real_, n_iter, 5L,
                         dimnames = list(NULL, c("b", "d", "g", "tau2", "sigma2")))
    accepted_flips <- 0L
    attempted_flips <- 0L
    accepted_domains <- 0L
    attempted_domains <- 0L

    for (it in seq_len(as.integer(n_iter))) {
      # STEP 1: Update d,g with variance integrated out (and beta also
      # integrated in default mode). The legacy no-design path preserves
      # its intermediate unused tau2 draw for RNG compatibility.
      if (is.null(fixed_gp)) {
        if (b == 1L) {
          gp <- bmim_gp_update(y, xi, D_by_loc, gp, gp_prior,
                               log_marginal = log_marginal,
                               draw_scale = is.null(scalar))
        } else {
          gp$d <- bmim_joint_sample_d_prior(gp_prior)
          gp$g <- rgamma(1L, shape = gp_prior$g_shape,
                        rate = gp_prior$g_rate)
        }
      }

      # STEP 2: Update a limited number of domains per sweep.
      # This number, not the number of candidate windows, controls how many
      # domain proposals are evaluated. Under inactive xi or b, sample prior.
      for (move in seq_len(as.integer(n_domain_updates_per_sweep))) {
        p <- sample.int(P, 1L)
        upd <- bmim_joint_update_domain(
          p, xi, A, D_by_loc, D_total, X_s, loc_slices, domain_priors,
          y, b, gp, gp_prior, domain_prior_jump_prob, log_alt = log_alt_kernel
        )
        A <- upd$A
        D_by_loc <- upd$D_by_loc
        D_total <- upd$D_total
        attempted_domains <- attempted_domains + 1L
        accepted_domains <- accepted_domains + as.integer(upd$accepted)
      }

      # STEP 3: Joint activation/deactivation and interval proposals.
      # Each node is visited once, in a random order to avoid scan bias.
      for (p in sample.int(P, P)) {
        move <- bmim_joint_flip_node(
          p, xi, A, D_by_loc, D_total, X_s, loc_slices, domain_priors,
          y, b, gp, gp_prior, W_mat, p_prop, log_alt = log_alt_kernel
        )
        xi <- move$xi
        A <- move$A
        D_by_loc <- move$D_by_loc
        D_total <- move$D_total
        attempted_flips <- attempted_flips + as.integer(move$attempted)
        accepted_flips <- accepted_flips + as.integer(move$accepted)
      }

      # Rebuild the selected distance once per sweep. This is inexpensive
      # relative to GP factorizations and prevents incremental subtraction
      # roundoff from accumulating across many accepted interval moves.
      D_total <- bmim_gp_selected_distance(D_by_loc, xi)

      # STEP 4: Optional global noise-only versus GP Bernoulli Gibbs draw.
      # The null uses its own sigma0^2; both variance scales are integrated.
      if (global_indicator) {
        log_alt <- if (is.null(scalar)) bmim_global_log_integrated_tau(
          y, bmim_gp_covariance(D_total, gp$d, gp$g), gp_prior
        ) else log_alt_kernel(y, D_total, gp, gp_prior)
        if (sampling_beta) {
          log_null <- bmim_scalar_log_evidence(
            bmim_scalar_null_terms(y, scalar, beta),
            gp_prior$null_sigma_shape, gp_prior$null_sigma_rate)
        }
        b <- rbinom(1L, 1L,
                    prob = plogis(qlogis(pi_b) + log_alt - log_null))
      }

      # STEP 5: Restore the two variance scales by exact conditional draws.
      scales <- bmim_joint_draw_scales(b, y, D_total, gp, gp_prior, scalar, beta)
      # Restore beta only AFTER the active variance draw. In integrated mode
      # optional draws are output-only and never fed back into the chain.
      if (recording_beta) {
        beta_draw <- bmim_scalar_draw_beta(y, b, D_total, gp, scales, scalar)
        beta_samples[it, ] <- beta_draw
        if (sampling_beta) beta <- beta_draw
      }
      gp$tau2 <- scales$tau2
      xi_samples[it, ] <- xi
      starts[it, ] <- A[, 1L]
      ends[it, ] <- A[, 2L]
      parameters[it, ] <- c(b, gp$d, gp$g, gp$tau2, scales$sigma2)
    }
    chains[[ch]] <- list(xi = xi_samples, start = starts, end = ends,
                         parameters = parameters, beta = beta_samples,
                         diagnostics = c(node_acceptance =
                           if (attempted_flips) accepted_flips / attempted_flips else NA_real_,
                           domain_acceptance =
                           if (attempted_domains) accepted_domains / attempted_domains else NA_real_))
  }

  kept <- seq.int(as.integer(burnin) + 1L, as.integer(n_iter))
  b_post <- unlist(lapply(chains, function(ch) ch$parameters[kept, "b"]))
  xi_post <- do.call(rbind, lapply(chains, function(ch) ch$xi[kept, , drop = FALSE]))
  starts_post <- do.call(rbind, lapply(chains, function(ch) ch$start[kept, , drop = FALSE]))
  ends_post <- do.call(rbind, lapply(chains, function(ch) ch$end[kept, , drop = FALSE]))
  included <- matrix(0, nrow = P, ncol = max(vapply(loc_slices, length, integer(1L))))
  conditional <- included
  for (p in seq_len(P)) {
    for (t in seq_len(length(loc_slices[[p]]))) {
      event <- xi_post[, p] == 1L & starts_post[, p] <= t & ends_post[, p] >= t
      included[p, t] <- mean(event & b_post == 1L)
      conditional[p, t] <- if (sum(b_post == 1L)) mean(event[b_post == 1L]) else NA_real_
    }
  }
  # Padding for shorter nodes is NA, not a fictitious zero posterior.
  for (p in seq_len(P)) {
    if (length(loc_slices[[p]]) < ncol(included)) {
      idx <- seq.int(length(loc_slices[[p]]) + 1L, ncol(included))
      included[p, idx] <- NA_real_
      conditional[p, idx] <- NA_real_
    }
  }
  structure(list(
    chains = chains, posterior_gp_included = mean(b_post),
    # Auxiliary xi remains defined under b=0. The unconditional xi summary
    # therefore includes null-model prior draws and is NOT evidence of an
    # active node. Use the conditional and joint summaries for interpretation.
    posterior_node_included = colMeans(xi_post),
    posterior_gp_and_node_included = colMeans(
      xi_post * matrix(b_post, nrow = nrow(xi_post), ncol = P)
    ),
    posterior_node_included_given_gp = if (any(b_post == 1L))
      colMeans(xi_post[b_post == 1L, , drop = FALSE]) else rep(NA_real_, P),
    posterior_node_time_joint_gp = included,
    posterior_node_time_given_gp = conditional,
    max_domain_length = as.integer(max_domain_length),
    min_domain_length = as.integer(min_domain_length),
    domain_priors = domain_priors,
    gp_prior = gp_prior,
    W_mat = W_mat, p_prop = p_prop, pi_b = pi_b,
    fixed_gp = fixed_gp,
    beta_mode = if (is.null(scalar)) "none" else beta_mode,
    beta_prior_mean = if (is.null(scalar)) NULL else scalar$mu,
    beta_prior_cov = if (is.null(scalar)) NULL else scalar$B,
    beta_design_names = if (is.null(scalar)) NULL else scalar$names,
    beta_draws_saved = recording_beta,
    n_chains = as.integer(n_chains), n_iter = as.integer(n_iter),
    burnin = as.integer(burnin),
    interpretation = paste(
      "Joint node/domain MCMC with a fixed empirical-Bayes selection prior;",
      "the global probability is NOT a selection-adjusted frequentist p-value."
    )
  ), class = "bmim_joint_domain_result")
}

# Convenience wrapper: preserve legacy no-Z preprocessing; with Z, fit
# scalar effects coherently in the likelihood, never by OLS residualization.
# p_prop must be EXPLICIT. For the flat-VB benchmark the caller may choose
# rep(0.5,P); for EB-VB use the fixed probabilities from the VB stage.
# This wrapper assumes each node's full-time columns are adjacent, as in ISGM.
bmim_joint_domains_from_features <- function(path_features, W_mat, y,
                                              grps_vec, p_prop,
                                              max_domain_length,
                                              Z = NULL, connec = 1,
                                              beta_mode = c("integrated", "sampled"),
                                              response_scale = NULL,
                                              ...) {
  X <- as.matrix(path_features)
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) stop("Nonfinite path_features.")
  lens <- vapply(grps_vec, function(z) length(unlist(z)), integer(1L))
  if (any(lens < 1L) || sum(lens) != ncol(X)) {
    stop("grps_vec must describe all full-time columns in node order.")
  }
  boundaries <- cumsum(lens)
  loc_slices <- lapply(seq_along(lens), function(p) {
    seq.int(boundaries[p] - lens[p] + 1L, boundaries[p])
  })
  X_s <- scale(X)
  beta_mode <- match.arg(beta_mode)
  if (!is.numeric(y) || length(y) != nrow(X) || any(!is.finite(y)))
    stop("y must be finite numeric with one value per subject.")
  if (any(!is.finite(X_s))) stop("Check constant functional columns.")
  # Preserve old production targets when no scalar design was supplied.
  # For scalar models use original response units by default, or a fixed
  # positive scale supplied independently of the observed response.
  if (!is.null(response_scale) && (!is.numeric(response_scale) ||
      length(response_scale) != 1L || !is.finite(response_scale) ||
      response_scale <= 0)) stop("response_scale must be positive and finite.")
  legacy_scaling <- is.null(Z) && is.null(response_scale)
  y_s <- if (legacy_scaling) as.numeric(scale(y)) else
    as.numeric(y)/(if (is.null(response_scale)) 1 else response_scale)
  if (any(!is.finite(y_s))) stop("Nonfinite response after scaling.")
  design <- NULL
  z_center <- z_scale <- NULL
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) || nrow(Z) != nrow(X) || any(!is.finite(Z)))
      stop("Z must contain finite numeric covariates, one row per subject.")
    if (ncol(Z)) {
      Z_s <- scale(Z)
      if (any(!is.finite(Z_s))) stop("Remove constant scalar columns; intercept is added automatically.")
      z_center <- attr(Z_s, "scaled:center")
      z_scale <- attr(Z_s, "scaled:scale")
      if (is.null(colnames(Z_s))) colnames(Z_s) <- paste0("Z", seq_len(ncol(Z_s)))
      design <- cbind("(Intercept)" = 1, Z_s)
    } else design <- matrix(1, nrow(X), 1L, dimnames=list(NULL,"(Intercept)"))
  }
  fit <- bmim_joint_domain_mcmc(
    X_s = X_s, y_tilde = y_s, loc_slices = loc_slices,
    W_mat = W_mat * connec, p_prop = p_prop,
    max_domain_length = max_domain_length, Z = design,
    beta_mode = beta_mode, ...
  )
  fit$response_preprocessing <- list(
    method = if (legacy_scaling) "legacy_standardized" else "fixed_scale",
    center = if (legacy_scaling) mean(y) else 0,
    scale = if (legacy_scaling) sd(y) else if (is.null(response_scale)) 1 else response_scale,
    scalar_center = z_center, scalar_scale = z_scale)
  fit
}

print.bmim_joint_domain_result <- function(x, ...) {
  cat("BMIM joint node/domain MCMC (conditional empirical Bayes)\n")
  cat("  Max domain length:", x$max_domain_length,
      "min domain length:", x$min_domain_length, "\n")
  cat("  Scalar effects:", x$beta_mode, "\n")
  cat("  Chains:", x$n_chains, "iterations/chain:", x$n_iter,
      "burn-in:", x$burnin, "\n")
  cat("  P(GP included | y):", signif(x$posterior_gp_included, 5), "\n")
  cat("  Node inclusion conditional on GP (b=1):",
      paste(signif(x$posterior_node_included_given_gp, 3),
            collapse = ", "), "\n")
  cat("  WARNING:", x$interpretation, "\n")
  invisible(x)
}
