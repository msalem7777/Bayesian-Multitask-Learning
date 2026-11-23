# ============================================================================
# bmim_global_gp_test.R -- conditional empirical-Bayes GP-vs-no-GP comparison
#
# An independent module for the existing BMIM implementation.
# Two proper, distinguishable variance models:
#   b=0: y|sigma0^2 ~ N(0,sigma0^2*I), sigma0^2~IG(a0,r0)
#   b=1: y|xi,d,g,tau2 ~ N(0,tau2*[K_xi(d)+g*I]), tau2~IG(a1,r1)
# The signal amplitude tau2 and relative nugget g appear only under b=1.
# Integrate sigma0^2 and tau2 analytically; integrate d,g for b=1 by
# common-prior Monte Carlo. Both variance priors are proper.
#
# Within a RESPONSE-SELECTED SHORTLIST S, set the alternative model weights
# proportional to exp(Ising(xi)) * product(VB Bernoulli probabilities).
# The weights are normalized WITHIN S. Selection of S, VB probabilities,
# and functional domains uses the same response: the reported result is a
# CONDITIONAL, DATA-ADAPTIVE EMPIRICAL-BAYES MODEL COMPARISON, not a calibrated
# global Bayes factor or frequentist p-value. The full alternative evidence
# is NOT recovered by silently normalizing within S. Do not use its posterior
# null probability as a significance p-value.
#
# Dependencies: base R and bmim_gp_updates.R, sourced by caller.
# ============================================================================

bmim_logsumexp <- function(x) {
  if (length(x) == 0L) stop("Cannot sum an empty vector.")
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

bmim_logmeanexp <- function(x) bmim_logsumexp(x) - log(length(x))

bmim_shortlist <- function(samples, M, threshold = 0.05,
                           mass_target = 0.99) {
  samples <- as.matrix(samples)
  if (!nrow(samples) || !ncol(samples) || anyNA(samples) ||
      any(!(samples %in% c(0, 1))) || any(rowSums(samples) == 0)) {
    stop("samples must be nonempty binary configurations with no empty models.")
  }
  if (length(M) != 1L || !is.finite(M) || M < 1 || M != as.integer(M) ||
      length(threshold) != 1L || !is.finite(threshold) ||
      threshold < 0 || threshold > 1 ||
      length(mass_target) != 1L || !is.finite(mass_target) ||
      mass_target <= 0 || mass_target > 1) {
    stop("M must be a positive integer; threshold in [0,1]; mass_target in (0,1].")
  }
  keys <- apply(samples, 1L, paste0, collapse = "")
  unique_keys <- unique(keys)
  ids <- match(keys, unique_keys)
  counts <- tabulate(ids, nbins = length(unique_keys))
  ranks <- order(-counts, unique_keys)
  ranked_keys <- unique_keys[ranks]
  ranked_prob <- counts[ranks] / nrow(samples)
  K <- min(length(ranked_keys), max(
    as.integer(M),
    sum(ranked_prob > threshold),
    which(cumsum(ranked_prob) + 10 * .Machine$double.eps >= mass_target)[1L]
  ))
  chosen <- ranked_keys[seq_len(K)]
  row_ids <- match(chosen, keys)
  configurations <- samples[row_ids, , drop = FALSE]
  rownames(configurations) <- chosen
  list(
    configurations = configurations,
    probabilities = setNames(ranked_prob[seq_len(K)], chosen),
    covered_mass = sum(ranked_prob[seq_len(K)]),
    omitted_mass = 1 - sum(ranked_prob[seq_len(K)]),
    all_probabilities = setNames(ranked_prob, ranked_keys),
    M = M, threshold = threshold, mass_target = mass_target
  )
}

bmim_global_log_model_prior <- function(configurations, W_mat, p_prop) {
  configurations <- as.matrix(configurations)
  P <- ncol(configurations)
  if (!all(dim(W_mat) == c(P, P)) || length(p_prop) != P ||
      any(!is.finite(p_prop)) || any(p_prop <= 0 | p_prop >= 1)) {
    stop("Incompatible network and empirical-Bayes prior dimensions/probabilities.")
  }
  logw <- vapply(seq_len(nrow(configurations)), function(i) {
    xi <- as.numeric(configurations[i, ])
    spin <- 2 * xi - 1
    0.5 * as.numeric(crossprod(spin, W_mat %*% spin)) +
      sum(dbinom(xi, 1, p_prop, log = TRUE))
  }, numeric(1))
  logw - bmim_logsumexp(logw)
}

# Exact tau2 integration, retaining all normalizing constants common to b=0,1.
# If y|tau2,C ~ N(0,tau2*C), tau2~IG(a,r), then
# log p(y|C)=lgamma(a+n/2)-lgamma(a)+a*log(r) -n/2*log(2*pi)
#             -1/2 log|C| -(a+n/2)*log(r+y'C^{-1}y/2).
bmim_global_log_integrated_tau <- function(y, covariance, prior) {
  y <- as.numeric(y)
  terms <- bmim_gp_cholesky_terms(covariance, y)
  if (is.null(terms)) return(-Inf)
  a <- prior$tau_shape
  r <- prior$tau_rate
  n <- length(y)
  lgamma(a + n / 2) - lgamma(a) + a * log(r) -
    (n / 2) * log(2 * pi) - terms$logdet / 2 -
    (a + n / 2) * log(r + terms$quad / 2)
}

# An identifiable null: sigma0^2 alone scales white noise, without g or
# GP signal variance tau2. IG(shape, rate) uses the inverse variance scale.
bmim_global_null_logevidence <- function(y, prior) {
  n <- length(y)
  a0 <- prior$null_sigma_shape
  r0 <- prior$null_sigma_rate
  if (is.null(a0) || is.null(r0) || a0 <= 0 || r0 <= 0) {
    stop("The null requires proper positive sigma0^2 inverse-gamma parameters.")
  }
  lgamma(a0 + n/2) - lgamma(a0) + a0*log(r0) -
    (n/2)*log(2*pi) - (a0+n/2)*log(r0 + sum(y*y)/2)
}

bmim_global_draw_prior <- function(prior, n_draws) {
  if (length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 2 || n_draws != as.integer(n_draws)) {
    stop("n_draws must be an integer >= 2.")
  }
  component <- sample.int(length(prior$d_mix_weights),
                          n_draws, replace = TRUE,
                          prob = prior$d_mix_weights)
  d <- prior$distance_scale * vapply(seq_len(n_draws), function(i) {
    rgamma(1L, shape = prior$d_mix_shapes[component[i]],
           rate = prior$d_mix_rates[component[i]])
  }, numeric(1))
  g <- rgamma(n_draws, shape = prior$g_shape, rate = prior$g_rate)
  list(d = d, g = g)
}

# A likelihood-concentration diagnostic for prior-MC evidence integration.
# This is not MCMC ESS and is not a proof that evidence is accurate.
bmim_global_importance_ess <- function(log_values) {
  weights <- exp(log_values - bmim_logsumexp(log_values))
  1 / sum(weights^2)
}

bmim_global_gp_test <- function(selection_samples, y_tilde, D_by_loc,
                                gp_prior, W_mat, p_prop,
                                M, threshold = 0.05,
                                mass_target = 0.99, pi_b,
                                n_prior_draws = 512L, seed = 20261122L,
                                prior_draws = NULL) {
  y <- as.numeric(y_tilde)
  if (length(y) != nrow(D_by_loc[[1L]]) ||
      any(!is.finite(y)) || any(!is.finite(W_mat))) {
    stop("Invalid response, distance matrices, or network.")
  }
  if (length(pi_b) != 1L || !is.finite(pi_b) ||
      pi_b <= 0 || pi_b >= 1) stop("pi_b must be between 0 and 1.")

  selected <- bmim_shortlist(selection_samples, M, threshold, mass_target)
  configs <- selected$configurations
  if (ncol(configs) != length(D_by_loc)) {
    stop("Configurations and distance matrices disagree.")
  }
  log_model_prior <- bmim_global_log_model_prior(configs, W_mat, p_prop)
  model_prior <- exp(log_model_prior)

  if (is.null(prior_draws)) {
    # Local seed without permanently changing the caller's RNG state.
    rng_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (rng_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (rng_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
    prior_draws <- bmim_global_draw_prior(gp_prior, as.integer(n_prior_draws))
  }
  d <- prior_draws$d
  g <- prior_draws$g
  if (length(d) < 2L || length(d) != length(g) ||
      any(!is.finite(d)) || any(!is.finite(g)) ||
      any(d <= 0) || any(g <= 0)) {
    stop("prior_draws$d and prior_draws$g must be positive matching vectors.")
  }
  n_draws <- length(d)
  n <- length(y)
  K <- nrow(configs)
  log_alternative_draws <- matrix(-Inf, nrow = n_draws, ncol = K)
  colnames(log_alternative_draws) <- rownames(configs)
  for (k in seq_len(K)) {
    D <- bmim_gp_selected_distance(D_by_loc, as.numeric(configs[k, ]))
    for (s in seq_len(n_draws)) {
      C <- bmim_gp_covariance(D, d[s], g[s])
      log_alternative_draws[s, k] <-
        bmim_global_log_integrated_tau(y, C, gp_prior)
    }
  }
  # The noise-only evidence is analytic and contains NO g or GP scale.
  log_evidence_null <- bmim_global_null_logevidence(y, gp_prior)

  log_marginal_configs <- apply(log_alternative_draws, 2L, bmim_logmeanexp)
  log_evidence_alt <- bmim_logsumexp(log_marginal_configs + log_model_prior)
  log_bayes_factor_10 <- log_evidence_alt - log_evidence_null
  posterior_b <- plogis(qlogis(pi_b) + log_bayes_factor_10)
  # Joint posterior probabilities of retained configurations conditional on
  # the GP being present, within the data-adaptive restricted model.
  conditional_post_xi <- exp(log_marginal_configs + log_model_prior -
                               log_evidence_alt)

  alt_mix_log_by_draw <- apply(sweep(log_alternative_draws, 2L,
                                     log_model_prior, "+"),
                                1L, bmim_logsumexp)
  # There is no Monte Carlo integration in the null evidence.
  null_ess <- Inf
  alt_ess <- bmim_global_importance_ess(alt_mix_log_by_draw)
  # A second diagnostic: compare log Bayes factors estimated from the
  # first and second halves of the SAME prior-draw sample.  The halves
  # are independent conditional on the chosen evidence-estimation setup.
  halves <- split(seq_len(n_draws),
                  rep(seq_len(2L), times = c(floor(n_draws / 2),
                                               n_draws - floor(n_draws / 2))))
  bf_by_half <- vapply(halves, function(idx) {
    log_alt_k <- apply(log_alternative_draws[idx, , drop = FALSE],
                       2L, bmim_logmeanexp)
    bmim_logsumexp(log_alt_k + log_model_prior) -
      log_evidence_null
  }, numeric(1))
  log_bf_half_difference <- abs(diff(unname(bf_by_half)))
  evidence_warning <- if (alt_ess < 0.1 * n_draws ||
                          log_bf_half_difference > 0.5) {
    "UNSTABLE ALTERNATIVE EVIDENCE: low prior-MC ESS and/or disagreement between halves; increase draws or use a stronger estimator."
  } else "No extreme prior-MC concentration or half-sample disagreement; convergence still needs checking."

  structure(list(
    shortlist = selected,
    configurations = configs,
    conditional_model_prior = setNames(model_prior, rownames(configs)),
    conditional_model_posterior_given_gp =
      setNames(conditional_post_xi, rownames(configs)),
    log_evidence_each_config = setNames(log_marginal_configs, rownames(configs)),
    log_evidence_null = log_evidence_null,
    log_evidence_alternative_truncated = log_evidence_alt,
    log_bayes_factor_10_truncated = log_bayes_factor_10,
    posterior_gp_included_truncated = posterior_b,
    posterior_noise_only_truncated = 1 - posterior_b,
    prior_integration_ess = c(null = null_ess, alternative = alt_ess),
    log_bf_by_half = bf_by_half,
    log_bf_half_difference = as.numeric(log_bf_half_difference),
    n_prior_draws = n_draws,
    prior_draws = prior_draws,
    evidence_warning = evidence_warning,
    interpretation = paste(
      "Conditional empirical-Bayes comparison on a response-selected shortlist;",
      "NOT a calibrated null probability or selection-integrated Bayes factor."
    )
  ), class = "bmim_global_gp_result")
}

# Convenience entrypoint for ISGM(..., return_details=TRUE).
bmim_global_gp_from_isgm <- function(details, burnin = 0L, ...) {
  if (!is.list(details) || is.null(details$samples) ||
      is.null(details$gp_prior) || is.null(details$p_prop)) {
    stop("Pass the result from ISGM(..., return_details=TRUE).")
  }
  samples <- as.matrix(details$samples)
  if (length(burnin) != 1L || !is.finite(burnin) ||
      burnin < 0 || burnin >= nrow(samples) ||
      burnin != as.integer(burnin)) {
    stop("burnin must be a valid number of initial sweeps to discard.")
  }
  keep <- seq.int(as.integer(burnin) + 1L, nrow(samples))
  bmim_global_gp_test(
    selection_samples = samples[keep, , drop = FALSE],
    y_tilde = details$y_tilde,
    D_by_loc = details$D_by_loc,
    gp_prior = details$gp_prior,
    W_mat = details$W_mat,
    p_prop = details$p_prop,
    ...
  )
}

print.bmim_global_gp_result <- function(x, ...) {
  cat("BMIM GP vs noise-only (conditional, adaptive shortlist)\n")
  cat("  Configurations retained:", nrow(x$configurations), "\n")
  cat("  Posterior frequency mass retained:",
      signif(x$shortlist$covered_mass, 4), "\n")
  cat("  log BF10 on shortlist:",
      signif(x$log_bayes_factor_10_truncated, 5), "\n")
  cat("  P(b=1 | data, selected shortlist):",
      signif(x$posterior_gp_included_truncated, 5), "\n")
  cat("  Prior-integration ESS:",
      paste(signif(x$prior_integration_ess, 4), collapse = ", "), "\n")
  cat("  Half-sample log-BF disagreement:",
      signif(x$log_bf_half_difference, 4), "\n")
  cat("  WARNING:", x$interpretation, "\n")
  cat("  Integration diagnostic:", x$evidence_warning, "\n")
  invisible(x)
}

# Obtain actual Bernoulli indicator draws (and conditional model-index draws)
# from the estimated *restricted conditional* posterior. No MCMC transition
# or new data is used here: these are direct categorical/Bernoulli samples
# from the model probabilities estimated by bmim_global_gp_test().
bmim_global_draw_indicator <- function(result, n = 1000L) {
  if (!inherits(result, "bmim_global_gp_result")) {
    stop("result must come from bmim_global_gp_test().")
  }
  if (length(n) != 1L || !is.finite(n) || n < 1 ||
      n != as.integer(n)) stop("n must be a positive integer.")
  b <- rbinom(n, size = 1L, prob = result$posterior_gp_included_truncated)
  config <- rep(NA_character_, n)
  active <- which(b == 1L)
  if (length(active)) {
    config[active] <- sample(
      names(result$conditional_model_posterior_given_gp),
      size = length(active), replace = TRUE,
      prob = result$conditional_model_posterior_given_gp
    )
  }
  data.frame(b = b, configuration = config)
}
