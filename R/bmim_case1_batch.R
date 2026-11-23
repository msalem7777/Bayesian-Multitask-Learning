# ============================================================================
# BMIM Case 1 Slurm array runner.
#
# This is an ADDITIVE execution driver; it does not modify the validated GP,
# Ising, global Bernoulli, or joint node/domain samplers.
#
# Execute from the repository's R/ directory:
#   Rscript --vanilla bmim_case1_batch.R --self-test
#   BMIM_TASK_ID=0 Rscript --vanilla bmim_case1_batch.R
#
# In production, Slurm provides SLURM_ARRAY_TASK_ID automatically. All other
# simulation controls are explicit BMIM_* environment variables, documented
# in hpc/README.md. Defaults below represent a PILOT, not publication settings.
# ============================================================================

bmim_batch_integer <- function(name, fallback, minimum = 0L) {
  raw <- Sys.getenv(name, unset = as.character(fallback))
  if (!grepl("^[0-9]+$", raw)) {
    stop(name, " must be an unsigned integer; received: ", raw)
  }
  value <- suppressWarnings(as.integer(raw))
  if (is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  value
}

bmim_batch_probability <- function(name, fallback) {
  raw <- Sys.getenv(name, unset = as.character(fallback))
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || !is.finite(value) ||
      value <= 0 || value >= 1) {
    stop(name, " must be a probability strictly between zero and one.")
  }
  value
}

bmim_batch_task_info <- function(task_id, per_scenario) {
  if (length(task_id) != 1L || is.na(task_id) || task_id < 0L ||
      length(per_scenario) != 1L || is.na(per_scenario) ||
      per_scenario < 1L || as.double(task_id) >= 2 * as.double(per_scenario)) {
    stop("Task index must be in 0,...,2*BMIM_REPLICATES_PER_SCENARIO-1.")
  }
  list(
    truth = as.integer(task_id %/% per_scenario),
    replicate = as.integer(task_id %% per_scenario + 1L)
  )
}

bmim_batch_validate_tag <- function(tag) {
  if (length(tag) != 1L || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", tag)) {
    stop("BMIM_RUN_TAG must contain only letters, digits, '.', '_', or '-'.")
  }
  tag
}

bmim_batch_selected_nodes <- function(p, seed, randomize, fixed_csv) {
  if (randomize == 1L) {
    # Separate, deterministic RNG stream for node identities. The adapter
    # reseeds the data generator itself, so this draw cannot shift its RNG.
    set.seed(as.integer(seed + 1000003L))
    return(sort(sample.int(p, size = 3L, replace = FALSE)))
  }
  values <- suppressWarnings(as.integer(strsplit(fixed_csv, ",", fixed = TRUE)[[1L]]))
  if (length(values) != 3L || anyNA(values) ||
      any(values < 1L | values > p) || length(unique(values)) != 3L) {
    stop("BMIM_CHOSEN_NODES must name exactly three distinct nodes in 1,...,P.")
  }
  values
}

bmim_batch_self_test <- function() {
  stopifnot(
    identical(bmim_batch_task_info(0L, 2L),
              list(truth = 0L, replicate = 1L)),
    identical(bmim_batch_task_info(1L, 2L),
              list(truth = 0L, replicate = 2L)),
    identical(bmim_batch_task_info(2L, 2L),
              list(truth = 1L, replicate = 1L)),
    identical(bmim_batch_task_info(3L, 2L),
              list(truth = 1L, replicate = 2L)),
    inherits(try(bmim_batch_task_info(4L, 2L), silent = TRUE), "try-error"),
    identical(bmim_batch_selected_nodes(10L, 19700L, 0L, "2,5,8"),
              c(2L, 5L, 8L)),
    identical(bmim_batch_selected_nodes(10L, 19700L, 1L, "2,5,8"),
              bmim_batch_selected_nodes(10L, 19700L, 1L, "2,5,8")),
    identical(bmim_batch_validate_tag("case1-pilot"), "case1-pilot")
  )
  cat("BMIM CASE 1 BATCH RUNNER SELF-TEST: PASS\n")
}

bmim_batch_main <- function() {
  # These are SAFETY defaults for the first benchmark, not finalized scientific
  # choices for the full paper. Override via explicit BMIM_* variables.
  n <- bmim_batch_integer("BMIM_N", 200L, 3L)
  p <- bmim_batch_integer("BMIM_P", 10L, 3L)
  max_length <- bmim_batch_integer("BMIM_LMAX", 20L, 1L)
  iterations <- bmim_batch_integer("BMIM_ITER", 1000L, 2L)
  burnin <- bmim_batch_integer("BMIM_BURNIN", 250L, 0L)
  chains <- bmim_batch_integer("BMIM_CHAINS", 4L, 1L)
  per_scenario <- bmim_batch_integer("BMIM_REPLICATES_PER_SCENARIO", 1L, 1L)
  seed_base <- bmim_batch_integer("BMIM_BASE_SEED", 19700L, 1L)
  pi_b <- bmim_batch_probability("BMIM_PI_B", 0.5)
  prior <- Sys.getenv("BMIM_PRIOR", unset = "vb")
  randomize <- bmim_batch_integer("BMIM_RANDOMIZE_NODES", 0L, 0L)
  tag <- bmim_batch_validate_tag(Sys.getenv("BMIM_RUN_TAG", unset = "pilot"))
  fixed_nodes <- Sys.getenv("BMIM_CHOSEN_NODES", unset = "2,5,8")
  task_id <- bmim_batch_integer("SLURM_ARRAY_TASK_ID",
                                bmim_batch_integer("BMIM_TASK_ID", 0L), 0L)

  if (iterations <= burnin) stop("BMIM_BURNIN must be less than BMIM_ITER.")
  if (!prior %in% c("vb", "flat")) stop("BMIM_PRIOR must be 'vb' or 'flat'.")
  if (!randomize %in% c(0L, 1L)) stop("BMIM_RANDOMIZE_NODES must be 0 or 1.")
  if (max_length > 100L) stop("This Case 1 generator has 100 time points per node.")
  if (as.double(seed_base) + 2 * per_scenario + 1000003 > .Machine$integer.max) {
    stop("Base seed plus task offsets exceeds R's integer range.")
  }
  info <- bmim_batch_task_info(task_id, per_scenario)
  seed <- as.integer(seed_base + task_id)
  chosen_nodes <- bmim_batch_selected_nodes(p, seed, randomize, fixed_nodes)

  output_base <- Sys.getenv("BMIM_RESULTS_DIR", unset =
                               file.path(Sys.getenv("HOME"), "BMIM_results"))
  if (!nzchar(output_base)) stop("BMIM_RESULTS_DIR must be nonempty.")
  run_dir <- file.path(
    output_base, tag,
    sprintf("case1_N%d_P%d_L%d", n, p, max_length),
    sprintf("iter%d_burn%d_ch%d_%s_pi%s_random%d",
            iterations, burnin, chains, prior,
            format(pi_b, scientific = FALSE, trim = TRUE), randomize)
  )
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(run_dir)) stop("Could not create run directory: ", run_dir)
  output_file <- file.path(run_dir,
      sprintf("truth%d_rep%04d_seed%d.rds", info$truth, info$replicate, seed))

  # Restart-safe: never recompute or silently overwrite a completed replicate.
  if (file.exists(output_file)) {
    previous <- tryCatch(readRDS(output_file), error = function(e) NULL)
    if (is.null(previous) || !is.list(previous) ||
        !identical(previous$batch_metadata$seed, seed) ||
        !identical(previous$batch_metadata$task_id, task_id)) {
      stop("Existing output is incomplete or does not match this task: ", output_file)
    }
    cat("SKIP: completed replicate already exists: ", output_file, "\n", sep = "")
    return(invisible(previous))
  }

  dependencies <- c("mvtnorm", "ClusterR")
  if (prior == "flat") dependencies <- "mvtnorm"
  missing <- dependencies[!vapply(dependencies, requireNamespace,
                                   logical(1L), quietly = TRUE)]
  if (length(missing)) stop("Missing R dependencies: ", paste(missing, collapse = ", "))

  # Source exactly the modules exercised in the passing integration smoke.
  if (!file.exists("bmim_joint_simulation.R")) {
    stop("Run this script from the repository's R/ directory.")
  }
  for (module in c("data_generation.R", "run_simulation.R",
                   "correlation_matrix_constructor.R", "ISGM_core.R",
                   "bmim_gp_updates.R", "bmim_global_gp_test.R",
                   "bmim_joint_domains.R", "bmim_joint_simulation.R")) {
    source(module, local = FALSE)
  }
  if (!exists("bmim_sim_valid_covariance", mode = "function")) {
    stop("Updated positive-definite data generator is required.")
  }

  cfg <- make_default_cfg(
    seed = NULL, sim_setting = "functional_linear", sim_iters = 1L,
    nx = n, num_locs = p, num_T = 20L, chosen_locs = chosen_nodes
  )
  metadata <- list(
    scenario = "functional_linear", truth = info$truth,
    replicate = info$replicate, task_id = task_id, seed = seed,
    chosen_nodes = chosen_nodes, n = n, p = p, max_domain_length = max_length,
    iterations = iterations, burnin = burnin, chains = chains,
    prior_mode = prior, pi_b = pi_b, randomize_nodes = randomize,
    job_id = Sys.getenv("SLURM_JOB_ID", unset = "not_slurm"),
    git_commit = tryCatch(system2("git", c("rev-parse", "HEAD"),
                                 stdout = TRUE, stderr = FALSE)[1L],
                          error = function(e) "unknown")
  )
  cat("START:", paste(vapply(names(metadata), function(key) {
    paste0(key, "=", paste(metadata[[key]], collapse = ","))
  }, character(1L)), collapse = " "), "\n")
  started <- Sys.time()
  result <- bmim_simulation_one(
    cfg = cfg, real_process = info$truth,
    max_domain_length = max_length, n_iter = iterations,
    burnin = burnin, n_chains = chains, prior_mode = prior,
    pi_b = pi_b, seed = seed
  )
  result$batch_metadata <- metadata
  result$batch_metadata$elapsed_seconds <- as.numeric(
    difftime(Sys.time(), started, units = "secs"))
  result$batch_metadata$decision_gp <-
    isTRUE(result$fit$posterior_gp_included > 0.5)

  # Save to a temporary sibling first; rename only after saveRDS succeeds.
  # A killed job therefore cannot leave a partially written final .rds file.
  temp_file <- paste0(output_file, ".tmp_", Sys.getpid())
  on.exit(if (file.exists(temp_file)) unlink(temp_file), add = TRUE)
  saveRDS(result, temp_file, compress = "gzip")
  if (file.exists(output_file)) stop("Output appeared during write: ", output_file)
  if (!file.rename(temp_file, output_file)) {
    stop("Failed to publish completed result atomically: ", output_file)
  }
  cat(sprintf(
    "DONE: truth=%d replicate=%d seed=%d P(GP)=%.6f decision(P>0.5)=%s elapsed=%.1f sec\n",
    info$truth, info$replicate, seed, result$fit$posterior_gp_included,
    if (result$batch_metadata$decision_gp) "GP" else "NULL",
    result$batch_metadata$elapsed_seconds
  ))
  cat("SAVED:", output_file, "\n")
  invisible(result)
}

if ("--self-test" %in% commandArgs(trailingOnly = TRUE)) {
  bmim_batch_self_test()
} else {
  # Turn unexpected warnings into errors rather than silently proceeding with
  # invalid synthetic covariances or numerical parameter estimates.
  options(warn = 2)
  bmim_batch_main()
}
