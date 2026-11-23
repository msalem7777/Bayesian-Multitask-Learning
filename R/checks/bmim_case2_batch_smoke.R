# Case-2 (historical R) CLI smoke, response-rule unit test, and short integration.
# Runs from R/: Rscript --vanilla -e "source('checks/bmim_case2_batch_smoke.R')"
# All numbers below are explicitly illustrative engineering test values.
# No publication results, original patient data, or figures are written.
script <- "bmim_case2_batch.R"
stopifnot(file.exists(script), file.exists("bmim_case2_threshold_response.R"))
output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                   c("--vanilla", script, "--self-test"),
                                   stdout = TRUE, stderr = TRUE))
status <- attr(output, "status")
if ((!is.null(status) && status != 0L) ||
    !any(grepl("BMIM CASE 2 (LEGACY) BATCH RUNNER SELF-TEST: PASS",
               output, fixed = TRUE))) {
  cat(paste(output, collapse = "\n"), "\n")
  stop("Case-2 CLI self-test failed.")
}
cat("CASE 2 BATCH CLI: PASS\n")

source("bmim_case2_threshold_response.R")
set.seed(913L)
K <- diag(c(1, -1, 2, -2))
alpha <- rep(1, 4)
rule <- bmim_case2_threshold_response(K, alpha,
                                      mu_positive = 5, sd_positive = 1,
                                      mu_negative = -5, sd_negative = 1)
stopifnot(identical(rule$group_positive, c(TRUE, FALSE, TRUE, FALSE)),
          identical(rule$mean, c(5, -5, 5, -5)),
          length(rule$y) == 4L, all(is.finite(rule$y)))
cat("MANUSCRIPT CASE-2 THRESHOLD RULE (GIVEN K, ALPHA): PASS\n")

for (module in c("data_generation.R", "run_simulation.R",
                 "correlation_matrix_constructor.R", "ISGM_core.R",
                 "bmim_gp_updates.R", "bmim_global_gp_test.R",
                 "bmim_joint_domains.R", "bmim_joint_simulation.R")) {
  source(module)
}
stats <- list(sim_mean_pat = 110, sim_sd_pat = 10,
              sim_mean_ctr = 100, sim_sd_ctr = 10)
stopifnot(identical(bmim_simulation_check_demo_stats(stats), stats),
          inherits(try(bmim_simulation_check_demo_stats(NULL), silent=TRUE),
                   "try-error"))
cfg <- make_default_cfg(seed = NULL, sim_setting = "real_or_noise",
                        sim_iters = 1L, nx = 24L, num_locs = 7L,
                        num_T = 20L, chosen_locs = c(2L, 5L, 7L))
for (truth in c(0L, 1L)) {
  fit <- bmim_simulation_one(
    cfg = cfg, real_process = truth, max_domain_length = 20L,
    n_iter = 40L, burnin = 10L, n_chains = 2L,
    prior_mode = "flat", pi_b = 0.5, seed = 91700L + truth,
    demo_stats = stats
  )
  stopifnot(identical(fit$config$sim_setting, "real_or_noise"),
            !isTRUE(fit$truth$true_domain_defined),
            length(fit$truth$true_timepoints) == 0L,
            all(is.finite(fit$fit$posterior_gp_included)),
            all(fit$fit$posterior_gp_included >= 0),
            all(fit$fit$posterior_gp_included <= 1),
            all(vapply(fit$fit$chains, function(ch) {
              all(ch$end - ch$start + 1L <= 20L) &&
                all(is.finite(ch$parameters))
            }, logical(1L))))
  cat("LEGACY CASE 2 INTEGRATION, truth=", truth, ": PASS\n", sep="")
}
cat("BMIM CASE 2 BATCH AND INTEGRATION SMOKE: PASS\n")
