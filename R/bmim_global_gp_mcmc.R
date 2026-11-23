# ============================================================================
# bmim_global_gp_mcmc.R -- conditional global GP-inclusion MCMC.
# Requires bmim_gp_updates.R and bmim_global_gp_test.R to be sourced first.
#
# CONDITIONAL model on the response-selected shortlist S:
#   b ~ Bernoulli(pi_b)
#   xi in S ~ normalized Ising(xi) * VB(xi), independent of b a priori
#   d ~ Gamma-mixture, g ~ Gamma(1,1), independent of b a priori
#   b=0: y | sigma0^2 ~ N(0, sigma0^2 I)
#        sigma0^2 ~ IG(null_sigma_shape,null_sigma_rate)
#   b=1: y | tau2,d,g,xi ~ N(0,tau2*(K_xi(d)+g I))
#        tau2 ~ IG(tau_shape,tau_rate)
#
# We integrate BOTH variance parameters out when updating (b,xi,d,g),
# then sample the active model's variance from its conjugate conditional.
# This prevents the unidentifiable tau2*g product under b=0. Under b=0,
# d,g,xi and tau2 are independent auxiliary prior variables; under b=1,
# sigma0^2 is independent auxiliary prior variable. These auxiliary draws
# ensure the jump back to either model uses a coherent joint prior.
#
# THIS IS A CONDITIONAL EMPIRICAL-BAYES COMPARISON: the shortlist, domains,
# and VB factor were constructed from the response. Do not interpret the
# posterior null frequency as a selection-adjusted frequentist p-value.
# ============================================================================

# Conditional white-noise Gaussian log likelihood (common constant omitted).
# This function uses sigma0^2 directly; there is no g or GP scale here.
bmim_global_null_loglik <- function(y, sigma2) {
  if (length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0) {
    return(-Inf)
  }
  n <- length(y)
  -0.5 * (n * log(sigma2) + sum(y*y) / sigma2)
}

# Draw the GP length scale independently from its two-component Gamma prior.
bmim_global_sample_d <- function(prior) {
  part <- sample.int(length(prior$d_mix_weights), 1L,
                     prob = prior$d_mix_weights)
  prior$distance_scale * rgamma(
    1L, shape = prior$d_mix_shapes[part], rate = prior$d_mix_rates[part]
  )
}

# Draw the GP nugget independently from its Gamma prior. Under the null,
# its value is auxiliary and cannot affect the observed-data likelihood.
bmim_global_sample_g <- function(prior) {
  rgamma(1L, shape = prior$g_shape, rate = prior$g_rate)
}

# Exact inverse-gamma draw. If v~IG(a,r), then 1/v~Gamma(a,rate=r).
bmim_global_inverse_gamma_draw <- function(shape, rate) {
  1 / rgamma(1L, shape = shape, rate = rate)
}

# Full conditional for the active model's scale; inactive variance is
# independently sampled from its prior for an explicit augmented state.
bmim_global_sample_scales <- function(b, y, D, gp, prior) {
  n <- length(y)
  if (b == 0L) {
    # Null variance is the only scale in the null likelihood.
    sigma2 <- bmim_global_inverse_gamma_draw(
      prior$null_sigma_shape + n / 2,
      prior$null_sigma_rate + sum(y*y) / 2
    )
    # Auxiliary GP scale is independent of y when the GP is absent.
    tau2 <- bmim_global_inverse_gamma_draw(
      prior$tau_shape, prior$tau_rate
    )
  } else {
    # GP scale conditional on K_xi(d)+g I.
    terms <- bmim_gp_cholesky_terms(
      bmim_gp_covariance(D, gp$d, gp$g), y
    )
    if (is.null(terms)) stop("GP covariance failed Cholesky decomposition.")
    tau2 <- bmim_global_inverse_gamma_draw(
      prior$tau_shape + n / 2,
      prior$tau_rate + terms$quad / 2
    )
    # Auxiliary null scale is independent of y when GP is present.
    sigma2 <- bmim_global_inverse_gamma_draw(
      prior$null_sigma_shape, prior$null_sigma_rate
    )
  }
  list(tau2 = tau2, sigma2 = sigma2)
}

bmim_global_gp_mcmc <- function(selection_samples, y_tilde, D_by_loc,
                                gp_prior, W_mat, p_prop, M, pi_b,
                                n_iter, burnin, n_chains,
                                threshold = 0.05, mass_target = 0.99,
                                seed = 20261122L) {
  if (length(n_iter) != 1L || !is.finite(n_iter) ||
      n_iter < 2L || n_iter != as.integer(n_iter) ||
      length(burnin) != 1L || !is.finite(burnin) ||
      burnin < 0L || burnin >= n_iter ||
      burnin != as.integer(burnin) ||
      length(n_chains) != 1L || !is.finite(n_chains) ||
      n_chains < 1L || n_chains != as.integer(n_chains)) {
    stop("Require integer n_iter>=2, 0<=burnin<n_iter, n_chains>=1.")
  }
  if (length(pi_b) != 1L || !is.finite(pi_b) ||
      pi_b <= 0 || pi_b >= 1) {
    stop("pi_b must be strictly inside (0,1).")
  }
  y <- as.numeric(y_tilde)
  if (length(y) != nrow(D_by_loc[[1L]]) || any(!is.finite(y))) {
    stop("Response and distance-matrix dimensions must agree.")
  }
  if (is.null(gp_prior$null_sigma_shape) ||
      is.null(gp_prior$null_sigma_rate)) {
    stop("Supply a gp_prior with proper independent null variance prior.")
  }
  selected <- bmim_shortlist(selection_samples, M = M,
                             threshold = threshold,
                             mass_target = mass_target)
  configs <- selected$configurations
  if (ncol(configs) != length(D_by_loc)) {
    stop("Configuration dimensions disagree.")
  }
  model_prior <- exp(bmim_global_log_model_prior(configs, W_mat, p_prop))
  K <- nrow(configs)
  D_list <- lapply(seq_len(K), function(k) {
    bmim_gp_selected_distance(D_by_loc, as.numeric(configs[k, ]))
  })
  log_null_evidence <- bmim_global_null_logevidence(y, gp_prior)

  rng_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (rng_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (rng_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  chains <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    b <- as.integer(ch %% 2L == 0L)
    idx <- sample.int(K, 1L, prob = model_prior)
    gp <- bmim_gp_initialize(gp_prior)
    chain <- matrix(
      NA_real_, nrow = n_iter, ncol = 6L,
      dimnames = list(NULL, c("b", "model_index", "d", "g", "tau2", "sigma2"))
    )

    for (it in seq_len(n_iter)) {
      # 1. Model-configuration update with both variance scales integrated.
      if (b == 0L) {
        # Under the null no configuration enters the data likelihood.
        idx <- sample.int(K, 1L, prob = model_prior)
      } else if (K > 1L) {
        proposed_idx <- sample.int(K, 1L, prob = model_prior)
        if (proposed_idx != idx) {
          ll_current <- bmim_gp_log_marginal_likelihood(
            y, D_list[[idx]], gp$d, gp$g, gp_prior
          )
          ll_proposed <- bmim_gp_log_marginal_likelihood(
            y, D_list[[proposed_idx]], gp$d, gp$g, gp_prior
          )
          # q(index*)=prior(index*) cancels the Ising/VB prior ratio.
          if (is.finite(ll_proposed) &&
              log(runif(1L)) < min(0, ll_proposed - ll_current)) {
            idx <- proposed_idx
          }
        }
      }

      # 2. Update only parameters identifiable in the current model.
      if (b == 1L) {
        # This step integrates out tau2 while updating d,g; its tau2 draw
        # will be refreshed after the Bernoulli update below.
        gp <- bmim_gp_update(
          y, as.numeric(configs[idx, ]), D_by_loc, gp, gp_prior
        )
      } else {
        # Auxiliary GP parameters are prior draws under the null.
        gp$d <- bmim_global_sample_d(gp_prior)
        gp$g <- bmim_global_sample_g(gp_prior)
        gp$tau2 <- bmim_global_inverse_gamma_draw(
          gp_prior$tau_shape, gp_prior$tau_rate
        )
      }

      # 3. Collapsed Bernoulli update: BOTH variance scales integrated out.
      log_alt <- bmim_global_log_integrated_tau(
        y, bmim_gp_covariance(D_list[[idx]], gp$d, gp$g), gp_prior
      )
      if (!is.finite(log_alt)) stop("Nonfinite alternative evidence.")
      prob_b1 <- plogis(qlogis(pi_b) + log_alt - log_null_evidence)
      if (!is.finite(prob_b1)) stop("Nonfinite Bernoulli conditional.")
      b <- rbinom(1L, 1L, prob = prob_b1)

      # 4. Refresh both scales conditional on the CURRENT indicator state.
      # Only the active scale informs the likelihood; the other is
      # auxiliary and comes from its proper prior.
      scales <- bmim_global_sample_scales(
        b, y, D_list[[idx]], gp, gp_prior
      )
      gp$tau2 <- scales$tau2
      chain[it, ] <- c(b, idx, gp$d, gp$g, scales$tau2, scales$sigma2)
    }
    chains[[ch]] <- chain
  }

  post <- do.call(rbind, lapply(chains, function(z) {
    z[seq.int(burnin + 1L, nrow(z)), , drop = FALSE]
  }))
  b_post <- post[, "b"]
  model_prob <- vapply(seq_len(K), function(k) {
    mean(b_post == 1L & post[, "model_index"] == k)
  }, numeric(1))
  names(model_prob) <- rownames(configs)
  structure(list(
    posterior_gp_included_truncated = mean(b_post),
    posterior_noise_only_truncated = mean(b_post == 0L),
    joint_posterior_gp_configurations = model_prob,
    chains = chains,
    samples_after_burnin = post,
    shortlist = selected,
    configs = configs,
    conditional_model_prior = setNames(model_prior, rownames(configs)),
    log_null_evidence = log_null_evidence,
    n_chains = n_chains, n_iter = n_iter, burnin = burnin,
    interpretation = paste(
      "Conditional response-adaptive shortlisted empirical-Bayes comparison;",
      "not a selection-integrated global posterior or frequentist p-value."
    )
  ), class = "bmim_global_gp_mcmc_result")
}

bmim_global_gp_mcmc_from_isgm <- function(details, burnin_selection = 0L, ...) {
  if (!is.list(details) || is.null(details$samples) ||
      is.null(details$gp_prior)) {
    stop("Pass result from ISGM(..., return_details=TRUE).")
  }
  if (length(burnin_selection) != 1L || !is.finite(burnin_selection) ||
      burnin_selection < 0L ||
      burnin_selection >= nrow(details$samples) ||
      burnin_selection != as.integer(burnin_selection)) {
    stop("Invalid burnin_selection.")
  }
  kept <- seq.int(burnin_selection + 1L, nrow(details$samples))
  bmim_global_gp_mcmc(
    selection_samples = details$samples[kept, , drop = FALSE],
    y_tilde = details$y_tilde, D_by_loc = details$D_by_loc,
    gp_prior = details$gp_prior, W_mat = details$W_mat,
    p_prop = details$p_prop, ...
  )
}

print.bmim_global_gp_mcmc_result <- function(x, ...) {
  cat("BMIM Bernoulli GP-indicator MCMC (conditional EB model)\n")
  cat("  Configurations retained:", nrow(x$configs),
      "; posterior-frequency mass covered:",
      signif(x$shortlist$covered_mass, 4), "\n")
  cat("  Chains:", x$n_chains, "sweeps/chain:", x$n_iter,
      "burn-in:", x$burnin, "\n")
  cat("  P(b=1 | data, restricted model):",
      signif(x$posterior_gp_included_truncated, 5), "\n")
  cat("  P(b=0 | data, restricted model):",
      signif(x$posterior_noise_only_truncated, 5), "\n")
  cat("  NOTE:", x$interpretation, "\n")
  invisible(x)
}
