# =============================================================================
# bmim_joint_simulation.R -- Reproducible synthetic-data integration for BMIM.
#
# This adapter uses the modular vb2 prior constructor; the data generator
# and joint GP likelihood remain unchanged.
#
# Required source order (from the project's R/ directory):
#   source("data_generation.R")
#   source("run_simulation.R")
#   source("correlation_matrix_constructor.R")
#   source("ISGM_core.R")                 # legacy helpers (if needed)
#   source("bmim_gp_updates.R")
#   source("bmim_global_gp_test.R")
#   source("bmim_joint_domains.R")
#   source("bmim_joint_simulation.R")
#
# Scientific controls must be passed explicitly by the caller. A smoke-test
# driver exists in checks/bmim_joint_simulation_smoke.R.
# =============================================================================

# Split each node's full, contiguous time columns in the concatenated design.
bmim_simulation_full_design <- function(full_dat2) {
  if (!is.list(full_dat2) || length(full_dat2) < 1L) {
    stop("full_dat2 must be a nonempty list, one subject-by-time matrix per node.")
  }
  blocks <- lapply(full_dat2, as.matrix)
  n <- nrow(blocks[[1L]])
  lens <- vapply(blocks, ncol, integer(1L))
  if (n < 3L || any(lens < 1L) ||
      any(vapply(blocks, nrow, integer(1L)) != n) ||
      any(!vapply(blocks, is.numeric, logical(1L))) ||
      any(!vapply(blocks, function(x) all(is.finite(x)), logical(1L)))) {
    stop("Each node must contain finite numeric data for the same subjects.")
  }
  X <- do.call(cbind, blocks)
  boundaries <- cumsum(lens)
  loc_slices <- lapply(seq_along(lens), function(p) {
    seq.int(boundaries[p] - lens[p] + 1L, boundaries[p])
  })
  list(X = X, loc_slices = loc_slices,
       grps_vec = lapply(lens, seq_len), n_subjects = n)
}

# Construct empirical weights once on complete trajectories. Flat is a
# separately named benchmark; no mixture or association calculation is needed.
bmim_simulation_selection_prior <- function(X, y, loc_slices, W_mat,
                                            mode, vb_summary = "auto",
                                            Z = NULL, details = FALSE,
                                            vb_seed = NULL) {
  if (length(mode) != 1L || !mode %in% c("flat", "vb"))
    stop("Specify prior mode explicitly: 'flat' or 'vb'.")
  if (identical(mode, "flat")) {
    out <- list(weights=rep(0.5,length(loc_slices)), version="flat", association=NULL)
  } else {
    out <- bmim_vb_prior(X,y,loc_slices,W_mat,Z,vb_summary,seed=vb_seed)
  }
  if (details) out else out$weights
}

# Validate the four ORIGINAL Case-2 group moments. These are computed from
# an external ABIDE dataset, never committed to Git. The function refuses
# missing or degenerate group scales, rather than inserting guessed values.
bmim_simulation_check_demo_stats <- function(demo_stats) {
  names_required <- c("sim_mean_pat", "sim_sd_pat",
                      "sim_mean_ctr", "sim_sd_ctr")
  if (!is.list(demo_stats) ||
      !all(names_required %in% names(demo_stats))) {
    stop("Case 2 needs demo_stats with four group mean/SD fields.")
  }
  out <- demo_stats[names_required]
  if (!all(vapply(out, function(x) is.numeric(x) && length(x) == 1L &&
                  is.finite(x), logical(1L))) ||
      out$sim_sd_pat <= 0 || out$sim_sd_ctr <= 0) {
    stop("Case 2 group means must be finite and SD values positive.")
  }
  out
}

# Fit one data-generating replicate, without FocusGP or the legacy score test.
#
# The legacy 'functional_linear' generator fixes signal support to timepoints
# 61:80 of a 5*num_T trajectory at the three chosen nodes.
# The legacy 'real_or_noise' generator uses an ABIDE-moment-based response
# mixture. It is the original R code's historical Case-2 driver, BUT IS NOT
# an implementation of the manuscript's K_gen %*% alpha threshold equation.
# For that branch, no true time interval is specified by the generator.
#
# When W_mat=NULL, connectivity is estimated from X alone using the existing
# build_W_matrix(mode='data') procedure. It does NOT use chosen_locs or y.
# For a fixed structural network, pass W_mat explicitly (symmetric P-by-P).
#
# prior_mode must be 'vb' (response-derived empirical Bayes) or
# 'flat' (a benchmark with all Bernoulli weights equal to 0.5).
#
# pi_b must be given explicitly: a number strictly between 0 and 1 to test
# GP versus independent noise, or NULL for a GP-always-present model.
#
# Each replicate's random seed is controlled here; the old simulator's
# cfg$seed is disabled in a copy so repeated replicates cannot silently
# regenerate the same dataset.
bmim_simulation_one <- function(cfg, real_process, max_domain_length,
                                n_iter, burnin, n_chains, prior_mode,
                                pi_b, seed, W_mat = NULL,
                                network_scale_full = 0.4,
                                connec = 1,
                                min_domain_length = 1L,
                                n_domain_updates_per_sweep = 1L,
                                domain_prior_jump_prob = 0.5,
                                demo_stats = NULL, vb_summary = "auto") {
  required <- c("simulate_one_dataset", "build_W_matrix",
                "bmim_joint_domains_from_features")
  missing <- required[!vapply(required, exists, logical(1L),
                               mode = "function", inherits = TRUE)]
  if (length(missing)) {
    stop("Source prerequisite modules first: ", paste(missing, collapse = ", "))
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("The legacy functional-data simulator requires mvtnorm.")
  }
  if (!(cfg$sim_setting %in% c("functional_linear", "real_or_noise"))) {
    stop("The adapter supports 'functional_linear' or 'real_or_noise' only.")
  }
  if (identical(cfg$sim_setting, "real_or_noise")) {
    # The original Case-2 code derives these from ABIDE demographics.
    # No fabricated defaults are permitted for a scientific simulation.
    demo_stats <- bmim_simulation_check_demo_stats(demo_stats)
  }
  if (length(cfg$chosen_locs) != 3L ||
      any(cfg$chosen_locs < 1L | cfg$chosen_locs > cfg$num_locs) ||
      length(unique(cfg$chosen_locs)) != 3L ||
      cfg$num_T < 16L || cfg$nx < 3L) {
    stop("Legacy functional generators require 3 distinct chosen nodes, num_T >= 16, nx >= 3.")
  }
  if (length(real_process) != 1L || !real_process %in% c(0, 1) ||
      length(seed) != 1L || !is.finite(seed) || seed < 1 ||
      seed != as.integer(seed)) {
    stop("real_process must be 0/1 and seed a positive integer.")
  }
  if (!is.null(W_mat) && (!is.matrix(W_mat) ||
      !identical(dim(W_mat), c(cfg$num_locs, cfg$num_locs)) ||
      any(!is.finite(W_mat)) || max(abs(W_mat - t(W_mat))) > 1e-9)) {
    stop("A supplied structural W_mat must be a finite symmetric P-by-P matrix.")
  }
  if (is.null(W_mat) && (length(network_scale_full) != 1L ||
      !is.finite(network_scale_full))) {
    stop("network_scale_full must be a finite scalar.")
  }

  set.seed(as.integer(seed))
  cfg_local <- cfg
  cfg_local$seed <- NULL
  # For functional_linear, group moments are unused; for real_or_noise,
  # require the four actual group moments from an external source.
  moments <- if (identical(cfg$sim_setting, "real_or_noise")) demo_stats else
    list(sim_mean_pat = NA_real_, sim_sd_pat = NA_real_,
         sim_mean_ctr = NA_real_, sim_sd_ctr = NA_real_)
  sim <- simulate_one_dataset(
    cfg = cfg_local,
    sim_mean_pat = moments$sim_mean_pat, sim_sd_pat = moments$sim_sd_pat,
    sim_mean_ctr = moments$sim_mean_ctr, sim_sd_ctr = moments$sim_sd_ctr,
    real_process = as.integer(real_process)
  )
  design <- bmim_simulation_full_design(sim$full_dat2)
  if (length(sim$sim_scores) != design$n_subjects ||
      any(!is.finite(sim$sim_scores))) {
    stop("Generated continuous response is invalid.")
  }
  if (is.null(W_mat)) {
    W_mat <- build_W_matrix(
      mode = "data", full_dat2 = sim$full_dat2,
      t_length = ncol(sim$full_dat2[[1L]]),
      scale_full = network_scale_full
    )$W_func_full
  }
  if (!is.matrix(W_mat) || any(!is.finite(W_mat)) ||
      max(abs(W_mat - t(W_mat))) > 1e-9) {
    stop("The empirical network is not a finite symmetric matrix.")
  }
  vb <- bmim_simulation_selection_prior(
    design$X, sim$sim_scores, design$loc_slices, W_mat, prior_mode,
    vb_summary = vb_summary, details = TRUE
  )
  p_prop <- vb$weights

  # The wrapper standardizes X and y once, exactly like the ISGM preprocessing.
  # No response-dependent FocusGP window selection is performed.
  fit <- bmim_joint_domains_from_features(
    path_features = design$X, W_mat = W_mat, y = sim$sim_scores,
    grps_vec = design$grps_vec, p_prop = p_prop,
    max_domain_length = max_domain_length, connec = connec,
    n_iter = n_iter, burnin = burnin, n_chains = n_chains,
    min_domain_length = min_domain_length, pi_b = pi_b,
    n_domain_updates_per_sweep = n_domain_updates_per_sweep,
    domain_prior_jump_prob = domain_prior_jump_prob,
    seed = as.integer(seed + 1L)
  )
  list(
    fit = fit,
    truth = list(real_process = as.integer(real_process),
                 true_nodes = as.integer(sim$true_locs),
                 # Real-or-noise has NO defined domain truth in legacy R code.
                 true_domain_defined = identical(cfg$sim_setting, "functional_linear") &&
                   real_process == 1L,
                 true_timepoints = if (identical(cfg$sim_setting, "functional_linear") &&
                   real_process == 1L) 61:80 else integer(0)),
    config = list(seed = as.integer(seed), nx = cfg$nx,
                  sim_setting = cfg$sim_setting,
                  demo_stats = if (identical(cfg$sim_setting, "real_or_noise"))
                    demo_stats else NULL,
                  num_locs = cfg$num_locs, num_T = cfg$num_T,
                  max_domain_length = max_domain_length,
                  prior_mode = prior_mode, vb_summary = if (is.null(vb$association)) "flat" else vb$association$method,
                  vb_version = vb$version, pi_b = pi_b,
                  n_iter = n_iter, burnin = burnin, n_chains = n_chains),
    summary = data.frame(
      real_process = as.integer(real_process),
      seed = as.integer(seed),
      posterior_gp = fit$posterior_gp_included,
      prior_mode = prior_mode,
      max_domain_length = max_domain_length
    ),
    p_prop = p_prop,
    vb = vb,
    W_mat = W_mat
  )
}
