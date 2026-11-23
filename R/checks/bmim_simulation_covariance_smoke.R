# BMIM synthetic covariance regression smoke test.
# Run from R/: Rscript -e "source('checks/bmim_simulation_covariance_smoke.R')"
# The test is intentionally independent of joint-domain inference code.
files <- Filter(Negate(is.null),
                lapply(sys.frames(), function(frame) frame$ofile))
if (!length(files)) stop("Source this test from its actual .R path.")
root <- dirname(dirname(normalizePath(tail(files, 1L)[[1L]])))
source(file.path(root, "data_generation.R"))
source(file.path(root, "run_simulation.R"))
rm(files, root)

old_warn <- getOption("warn")
options(warn = 2)   # A warning from mvtnorm should fail this smoke test.

# These pairwise distances intentionally violate Euclidean geometry.
# The uncorrected correlation matrix must be indefinite.
D_bad <- matrix(c(0, 0.01, 80,
                  0.01, 0, 0.01,
                  80, 0.01, 0), nrow = 3L, byrow = TRUE)
raw <- 4 * exp(-0.5 * (D_bad - mean(D_bad)) / sd(as.vector(D_bad)))
raw_sd <- sqrt(diag(raw))
raw_cor <- raw / outer(raw_sd, raw_sd)
stopifnot(min(eigen(raw_cor, symmetric = TRUE,
                    only.values = TRUE)$values) < 0)

fixed <- bmim_sim_valid_covariance(D_bad, variance = 4)
stopifnot(isTRUE(all.equal(diag(fixed), diag(raw), tolerance = 1e-12)),
          isTRUE(all.equal(fixed, t(fixed), tolerance = 1e-12)),
          min(eigen(fixed, symmetric = TRUE,
                    only.values = TRUE)$values) > 0)
# One multiplier must be applied to all off-diagonal correlations.
idx <- upper.tri(raw)
ratios <- fixed[idx] / raw[idx]
stopifnot(max(ratios) - min(ratios) < 1e-10)

# When the original Gaussian covariance is valid, do not change it.
D_good <- as.matrix(dist(c(0, 1, 2)))^2
raw_good <- 4 * exp(-0.5 *
  (D_good - mean(D_good)) / sd(as.vector(D_good)))
fixed_good <- bmim_sim_valid_covariance(D_good, variance = 4)
stopifnot(isTRUE(all.equal(raw_good, fixed_good, tolerance = 1e-12)))
cat("COVARIANCE VALIDITY AND MINIMAL REPAIR: PASS\n")

# Exercise the real simulation generator with warnings promoted to errors.
cfg <- make_default_cfg(
  seed = NULL, sim_setting = "functional_linear", sim_iters = 1L,
  nx = 24L, num_locs = 7L, num_T = 20L,
  chosen_locs = c(2L, 5L, 7L)
)
for (truth in c(0L, 1L)) {
  set.seed(19500L + truth)
  sim <- simulate_one_dataset(
    cfg, NA_real_, NA_real_, NA_real_, NA_real_, real_process = truth
  )
  stopifnot(length(sim$full_dat2) == cfg$num_locs,
            all(vapply(sim$full_dat2, function(x)
              nrow(x) == cfg$nx && ncol(x) == cfg$t_length &&
              all(is.finite(as.matrix(x))), logical(1L))),
            all(is.finite(sim$sim_scores)))
}
options(warn = old_warn)
cat("SIMULATION NULL AND SIGNAL GENERATION: PASS\n")
cat("BMIM SIMULATION COVARIANCE SMOKE: PASS\n")
