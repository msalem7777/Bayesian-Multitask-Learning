# =============================================================================
# End-to-end SIMULATION integration smoke. This is not an accuracy study.
# It uses the legacy functional-linear generator but NEVER calls FocusGP.
# Run from R/: Rscript -e "source('checks/bmim_joint_simulation_smoke.R')"
# =============================================================================
files <- Filter(Negate(is.null),
                lapply(sys.frames(), function(frame) frame$ofile))
if (!length(files)) stop("Source this smoke test from its actual .R path.")
root <- dirname(dirname(normalizePath(tail(files, 1L)[[1L]])))
source(file.path(root, "data_generation.R"))
source(file.path(root, "run_simulation.R"))
source(file.path(root, "correlation_matrix_constructor.R"))
source(file.path(root, "ISGM_core.R"))
source(file.path(root, "bmim_gp_updates.R"))
source(file.path(root, "bmim_global_gp_test.R"))
source(file.path(root, "bmim_joint_domains.R"))
source(file.path(root, "bmim_joint_simulation.R"))
rm(files, root)

# Explicit engineering-only toy settings. NOT publication simulation settings.
cfg <- make_default_cfg(
  seed = NULL,
  sim_setting = "functional_linear",
  sim_iters = 1L,
  nx = 24L,
  num_locs = 7L,
  num_T = 20L,
  chosen_locs = c(2L, 5L, 7L)
)

# Verify the simulator/adapter preserves every node's FULL time course.
set.seed(19421)
standalone <- simulate_one_dataset(cfg, NA_real_, NA_real_, NA_real_,
                                   NA_real_, real_process = 1L)
design <- bmim_simulation_full_design(standalone$full_dat2)
stopifnot(nrow(design$X) == cfg$nx,
          ncol(design$X) == cfg$num_locs * cfg$t_length,
          length(design$loc_slices) == cfg$num_locs,
          all(vapply(design$loc_slices, length, integer(1L)) == cfg$t_length),
          identical(sort(as.integer(unlist(design$loc_slices))),
                    seq_len(ncol(design$X))))
cat("SIMULATION FULL-TRAJECTORY SHAPE: PASS\n")

# Test the NULL and SIGNAL branches, with flat Bernoulli weights explicitly
# identified as a baseline (not response-derived VB).
for (truth in c(0L, 1L)) {
  result <- bmim_simulation_one(
    cfg = cfg, real_process = truth,
    max_domain_length = 20L,
    n_iter = 40L, burnin = 10L, n_chains = 2L,
    prior_mode = "flat", pi_b = 0.5,
    seed = 19421L + truth,
    network_scale_full = 0.4,
    n_domain_updates_per_sweep = 1L
  )
  stopifnot(result$truth$real_process == truth,
            identical(result$truth$true_nodes,
                      as.integer(if (truth == 0L) rep(0, cfg$num_locs)
                                 else seq_len(cfg$num_locs) %in% cfg$chosen_locs)),
            length(result$p_prop) == cfg$num_locs,
            all(result$p_prop == 0.5),
            is.finite(result$fit$posterior_gp_included),
            result$fit$posterior_gp_included >= 0,
            result$fit$posterior_gp_included <= 1,
            all(is.finite(result$W_mat)),
            nrow(result$fit$posterior_node_time_joint_gp) == cfg$num_locs)
  for (chain in result$fit$chains) {
    stopifnot(all(chain$end >= chain$start),
              all(chain$end - chain$start + 1L <= 20L),
              all(rowSums(chain$xi) > 0L),
              all(is.finite(chain$parameters)))
  }
  cat(sprintf("SIMULATION %s INTEGRATION: PASS; posterior P(GP)=%.4f\n",
              if (truth == 0L) "NULL" else "SIGNAL",
              result$fit$posterior_gp_included))
}
cat("BMIM JOINT SIMULATION INTEGRATION SMOKE: PASS\n")
