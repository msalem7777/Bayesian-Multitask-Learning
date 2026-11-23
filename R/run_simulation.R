# =========================================================
# run_simulation.R
# End-to-end simulation driver (1 replicate or Monte Carlo)
# =========================================================
# This script is meant to be the *only* file you run directly.
# Everything else should live in R/ as package functions.
#
# Depends on (from earlier scripts):
#   - fn_dat_creator(), fn_dat_appendor(), fn_dat_preppendor()   [data generation]
#   - FocusGP(), FCV(), EI.search(), ...                         [domain selection]
#   - ISGM()                                                     [location selection]
#   - tau.score.test() (your original), plus estimate_rho(), etc  [tau tests]
#   - normalize() helper
#
# NOTE: I keep file IO optional. For an R package, prefer returning
# objects and letting the caller decide where/what to save.


# -------------------------------
# Configuration
# -------------------------------

#' Create a default config list
#' @export
make_default_cfg <- function(
  seed = NULL,
  sim_setting = c("real_or_noise", "functional_linear"),
  sim_iters = 100,
  nx = 200,
  num_locs = 10,
  num_T = 20,
  chosen_locs = c(2,5,7),
  chosen_one = 1,
  unchosen_one = 2,
  gt_sigma_u = 4,
  gt_sigma_pop = 0.25,
  gt_tau = 4,
  connec = 1.0,
  isgm_N = 300,
  focusgp_nfold = 5,
  focusgp_lambda = 0.05,
  focusgp_ninit = 12,
  focusgp_end = 30,
  lengthscale = NULL        # ✅ ADD THIS
){

  sim_setting <- match.arg(sim_setting)

  list(
    seed = seed,
    sim_setting = sim_setting,
    sim_iters = sim_iters,
    nx = nx,
    num_locs = num_locs,
    num_T = num_T,
    t_length = num_T * 5,

    chosen_locs = chosen_locs,
    chosen_one = chosen_one,
    unchosen_one = unchosen_one,

    gt_sigma_u = gt_sigma_u,
    gt_sigma_pop = gt_sigma_pop,
    gt_tau = gt_tau,

    connec = connec,
    isgm_N = isgm_N,

    focusgp_nfold = focusgp_nfold,
    focusgp_lambda = focusgp_lambda,
    focusgp_ninit = focusgp_ninit,
    focusgp_end = focusgp_end,
    lengthscale = lengthscale
  )
}


# -------------------------------
# Small helpers
# -------------------------------

#' Safe scale to matrix
#' @export
as_scaled <- function(x){
  scale(as.matrix(x))
}

#' Location prior adjustment (same logic as your script)
#' @export
prior_adjustment <- function(num_locs, nx){
  ifelse((num_locs / nx > 1), 2 * 3 * (num_locs / nx), 1)
}

#' Build a distance/affinity matrix with stronger ties among chosen_locs
#' (matches your generation logic structure)
#' @export
build_dist_mat <- function(num_locs, chosen_locs){
  dist_mat <- matrix(0, nrow = num_locs, ncol = num_locs)
  for(a in seq_len(num_locs)){
    for(b in seq_len(num_locs)){
      if(a < b){
        if((a %in% chosen_locs) && (b %in% chosen_locs)){
          dist_mat[a,b] <- dist_mat[b,a] <- rgamma(1, 4.5, 3.5)
        } else {
          dist_mat[a,b] <- dist_mat[b,a] <- rnorm(1, 4.5, 1.5)
        }
      }
    }
  }
  dist_mat
}

#' Build W matrices (Ising connectivity) like your script
#' @export
build_W_func <- function(num_locs, chosen_locs, t_length){
  W_func <- matrix(0, nrow = num_locs, ncol = num_locs)
  W_func_full <- matrix(0, nrow = num_locs, ncol = num_locs)

  for(l in seq_len(num_locs)){
    for(k in seq_len(num_locs)){
      if(l >= k) next
      if((l %in% chosen_locs) && (k %in% chosen_locs)){
        W_func[k,l] <- rnorm(1, 0.75, 0.05)
        W_func_full[k,l] <- 0.4 * W_func[k,l]
      } else {
        tmp <- rt(1, t_length - 2, 0)
        tmp <- tmp * ((1 - tmp^2) / (t_length - 2))
        W_func[k,l] <- tmp
        W_func_full[k,l] <- tmp
      }
      W_func[l,k] <- W_func[k,l]
      W_func_full[l,k] <- W_func_full[k,l]
    }
  }

  list(W_func = W_func, W_func_full = W_func_full)
}


# -------------------------------
# Data generation (single replicate)
# -------------------------------

#' Generate one simulated dataset (functional predictors by location)
#'
#' Returns:
#'  - full_dat2: list length num_locs, each a (nx x t_length) data.frame
#'  - sim_scores: numeric(nx)
#'  - true_status: 0/1 class labels per subject
#'  - true_locs: 0/1 relevance per location (only meaningful if real_process==1)
#'  - real_process: 0/1 indicator
#'  - dist_mat, chosen_locs
#'
#' @export
simulate_one_dataset <- function(cfg, sim_mean_pat, sim_sd_pat, sim_mean_ctr, sim_sd_ctr, real_process = NULL){

  if(!is.null(cfg$seed)) set.seed(cfg$seed)

  N <- cfg$nx
  num_locs <- cfg$num_locs
  num_T <- cfg$num_T
  t_length <- cfg$t_length

  chosen_locs <- cfg$chosen_locs
  chosen_one <- cfg$chosen_one
  unchosen_one <- cfg$unchosen_one

  # Decide whether this replicate has a "real process"
  if(is.null(real_process)){
    real_process <- rbinom(1, 1, 1.0)
  }

  true_locs <- rep(0, num_locs)
  if(real_process == 1) true_locs[chosen_locs] <- 1

  dist_mat <- build_dist_mat(num_locs, chosen_locs)

  # Create base functional data for 2 prototypes
  fn_dat_dfs <- fn_dat_creator(num_locs, 2, num_T, dist_mat = dist_mat, gt_tau = cfg$gt_tau)

  # Expand to N subjects
  for(i in 3:N){
    fn_dat_dfs <- append(fn_dat_dfs, list(matrix(0.0, nrow = num_locs, ncol = num_T)))
  }

  # Case/control indicator per subject
  true_status <- numeric(N)
  true_status[1:2] <- c(1, 0)

  for(i in 3:N){
    true_status[i] <- rbinom(1, 1, 0.5)

    if(real_process == 0){
      for(j in seq_len(num_locs)){
        fn_dat_dfs[[i]][j,] <- fn_dat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    } else {
      if(true_status[i] == 1){
        for(j in seq_len(num_locs)){
          if(j %in% chosen_locs){
            fn_dat_dfs[[i]][j,] <- fn_dat_dfs[[chosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
          } else {
            fn_dat_dfs[[i]][j,] <- fn_dat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
          }
        }
      } else {
        for(j in seq_len(num_locs)){
          fn_dat_dfs[[i]][j,] <- fn_dat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
        }
      }
    }

    # Make chosen_one nontrivial outside chosen_locs
    if(real_process == 0){
      for(j in seq_len(num_locs)){
        fn_dat_dfs[[chosen_one]][j,] <- fn_dat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    } else {
      for(j in seq_len(num_locs)){
        if(j %in% chosen_locs) next
        fn_dat_dfs[[chosen_one]][j,] <- fn_dat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    }
  }

  # Append / prepend blocks
  appendat_dfs <- fn_dat_appendor(num_locs, 2, num_T, fn_dat_dfs, dist_mat = dist_mat, gt_sigma_u = cfg$gt_sigma_u)
  for(i in 3:N) appendat_dfs <- append(appendat_dfs, list(matrix(0.0, nrow = num_locs, ncol = num_T)))

  for(i in seq_len(N)){
    if(i == unchosen_one) next
    for(j in seq_len(num_locs)){
      if((real_process == 1) && (j %in% chosen_locs) && (true_status[i] == 1)){
        appendat_dfs[[i]][j,] <- (runif(1, 0.950, 0.965)^(seq_len(num_T))) *
          (fn_dat_dfs[[i]][j, num_T/2] - appendat_dfs[[unchosen_one]][j, 1]) +
          appendat_dfs[[unchosen_one]][j,] +
          mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      } else {
        appendat_dfs[[i]][j,] <- appendat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    }
  }

  preppendat_dfs <- fn_dat_preppendor(num_locs, 2, num_T, fn_dat_dfs, dist_mat = dist_mat, gt_sigma_u = cfg$gt_sigma_u)
  for(i in 3:N) preppendat_dfs <- append(preppendat_dfs, list(matrix(0.0, nrow = num_locs, ncol = num_T)))
  for(i in seq_len(N)){
    if(i == unchosen_one) next
    for(j in seq_len(num_locs)){
      preppendat_dfs[[i]][j,] <- preppendat_dfs[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
    }
  }

  # Second cycle
  fn_dat_dfs2 <- fn_dat_creator(num_locs, 2, num_T, dist_mat = dist_mat, gt_tau = cfg$gt_tau)
  for(i in 3:N) fn_dat_dfs2 <- append(fn_dat_dfs2, list(matrix(0.0, nrow = num_locs, ncol = num_T)))

  for(i in 3:N){
    if(real_process == 0){
      for(j in seq_len(num_locs)){
        fn_dat_dfs2[[i]][j,] <- fn_dat_dfs2[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    } else {
      if(true_status[i] == 1){
        for(j in seq_len(num_locs)){
          if(j %in% chosen_locs){
            fn_dat_dfs2[[i]][j,] <- fn_dat_dfs2[[chosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
          } else {
            fn_dat_dfs2[[i]][j,] <- fn_dat_dfs2[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
          }
        }
      } else {
        for(j in seq_len(num_locs)){
          fn_dat_dfs2[[i]][j,] <- fn_dat_dfs2[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
        }
      }
    }

    if(real_process == 0){
      for(j in seq_len(num_locs)){
        fn_dat_dfs2[[chosen_one]][j,] <- fn_dat_dfs2[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    } else {
      for(j in seq_len(num_locs)){
        if(j %in% chosen_locs) next
        fn_dat_dfs2[[chosen_one]][j,] <- fn_dat_dfs2[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    }
  }

  appendat_dfs2 <- fn_dat_appendor(num_locs, 2, num_T, fn_dat_dfs2, dist_mat = dist_mat, gt_sigma_u = cfg$gt_sigma_u)
  for(i in 3:N) appendat_dfs2 <- append(appendat_dfs2, list(matrix(0.0, nrow = num_locs, ncol = num_T)))

  for(i in seq_len(N)){
    if(i == unchosen_one) next
    for(j in seq_len(num_locs)){
      if((real_process == 1) && (j %in% chosen_locs) && (true_status[i] == 1)){
        appendat_dfs2[[i]][j,] <- (runif(1, 0.950, 0.965)^(seq_len(num_T))) *
          (fn_dat_dfs2[[i]][j, num_T/2] - appendat_dfs2[[unchosen_one]][j, 1]) +
          appendat_dfs2[[unchosen_one]][j,] +
          mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      } else {
        appendat_dfs2[[i]][j,] <- appendat_dfs2[[unchosen_one]][j,] + mvtnorm::rmvnorm(1, mean = rep(0, num_T), sigma = cfg$gt_sigma_pop * diag(num_T))
      }
    }
  }

  # Combine blocks: (prepend) + (cycle1 + append1) + (cycle2 + append2)
  full_dat <- Map(cbind, fn_dat_dfs, appendat_dfs)
  full_dat <- Map(cbind, full_dat, fn_dat_dfs2)
  full_dat <- Map(cbind, full_dat, appendat_dfs2)
  full_dat <- Map(cbind, preppendat_dfs, full_dat)

  # Convert to per-location subject x time (nx x t_length) with random intercepts
  full_dat2 <- vector("list", length = num_locs)

  for(l in seq_len(num_locs)){
    full_dat2[[l]] <- as.data.frame(matrix(0, nrow = N, ncol = t_length))
  }

  for(l in seq_len(num_locs)){
    if((l %in% chosen_locs) && (real_process == 1)){
      rand_mu_1 <- sample(c(-1, 1), 1)
      rand_mu_0 <- -rand_mu_1
    } else {
      rand_mu_1 <- 0
      rand_mu_0 <- 0
    }

    for(i in seq_len(N)){
      full_dat2[[l]][i,] <- full_dat[[i]][l,] + rep(rnorm(1, (rand_mu_0 * (true_status[i] == 0)) + (rand_mu_1 * (true_status[i] == 1)), 1), t_length)
    }
  }

  # -------------------------------
  # Response generation (two simulation settings)
  # -------------------------------

  N_pat <- sum(true_status == 1)
  N_ctr <- sum(true_status == 0)

  if(cfg$sim_setting == "real_or_noise"){
    # Original mixture-based setting
    sim_scores <- (real_process == 1) * (
      rnorm(N, sim_mean_pat, sim_sd_pat) * (true_status == 1) +
        rnorm(N, sim_mean_ctr, sim_sd_ctr) * (true_status == 0)
    ) + (real_process == 0) * rnorm(N, sim_mean_ctr, sim_sd_ctr)

  } else if(cfg$sim_setting == "functional_linear"){
    # Functional linear model using true locations
    beta_1 <- beta_2 <- beta_3 <- rep(0, t_length)
    grid <- seq(0, 1, length.out = t_length)

    beta_1[61:80] <-  0.3/3 * dnorm(grid, mean = .675, sd = .03)[61:80]
    beta_2[61:80] <- -0.25/3 * dnorm(grid, mean = .70,  sd = .05)[61:80]
    beta_3[61:80] <- -0.40/3 * dnorm(grid, mean = .725, sd = .03)[61:80]

    if(real_process == 1){
      sim_scores <- 100 +
        as.matrix(full_dat2[[chosen_locs[1]]]) %*% beta_1 +
        as.matrix(full_dat2[[chosen_locs[2]]]) %*% beta_2 +
        as.matrix(full_dat2[[chosen_locs[3]]]) %*% beta_3 +
        rnorm(N, 0, 1)
    } else {
      sim_scores <- 100 + rnorm(N, 0, 1)
    }

    sim_scores <- as.numeric(sim_scores)
  }

  list(

    full_dat2 = full_dat2,
    sim_scores = sim_scores,
    true_status = true_status,
    true_locs = true_locs,
    real_process = real_process,
    chosen_locs = chosen_locs,
    dist_mat = dist_mat
  )
}


# -------------------------------
# Domain selection (FocusGP per location)
# -------------------------------

#' Select contiguous time domain per location using FocusGP
#'
#' Returns:
#'  - dat_trunc: list of truncated design matrices per location
#'  - grps_vec: list of selected column indices grouped by consecutive runs
#'  - domains: list of 0/1 length t_length per location (selected window)
#' @export
select_domains_focusgp <- function(full_dat2, sim_scores, cfg){

  num_locs <- cfg$num_locs
  t_length <- cfg$t_length
  N <- cfg$nx

  dat_trunc <- full_dat2
  grps_vec <- vector("list", length = num_locs)
  domains <- vector("list", length = num_locs)

  for(l in seq_len(num_locs)){

    X <- full_dat2[[l]][, 1:t_length, drop = FALSE]
    Z <- sim_scores

    # FocusGP can error because of GP fitting; let caller decide retry policy.
    GP_select_final <- FocusGP(
      X = X,
      Z = Z,
      nfold = cfg$focusgp_nfold,
      lambda = cfg$focusgp_lambda,
      ninit = cfg$focusgp_ninit,
      end = cfg$focusgp_end
    )

    a_sel <- GP_select_final$X[which.min(GP_select_final$maxei), 1]
    b_sel <- a_sel + GP_select_final$X[which.min(GP_select_final$maxei), 2]
    a_sel <- max(1, floor(a_sel))
    b_sel <- min(t_length, floor(b_sel))

    # Center each subject trajectory (like your scale(..., scale=FALSE))
    appnder <- as.data.frame(t(apply(dat_trunc[[l]][, 1:t_length, drop = FALSE], 1, scale, scale = FALSE)))
    colnames(appnder) <- paste0("t", seq_len(t_length))

    sel_idx <- a_sel:b_sel
    dat_trunc[[l]] <- appnder[, sel_idx, drop = FALSE]

    dom <- rep(0, t_length)
    dom[sel_idx] <- 1
    domains[[l]] <- dom

    consec <- split(sel_idx, cumsum(c(1, diff(sel_idx) != 1)))
    grps_vec[[l]] <- consec

    # optional intercept feature (kept separate if you want it)
    intrcpts <- rowMeans(full_dat2[[l]][, 1:t_length, drop = FALSE], na.rm = TRUE)
    dat_trunc[[l]]$intrcpts <- intrcpts
  }

  list(dat_trunc = dat_trunc, grps_vec = grps_vec, domains = domains)
}


# -------------------------------
# Build design matrices
# -------------------------------

#' Build numeric-safe integrated design matrix
#' @export
make_integrated_df <- function(dat_list){
  as.data.frame(do.call(cbind, dat_list))
}

#' Build full grouping vector in ISGM-compatible format
#' @export
make_full_grps_vec <- function(dat_list){
  lapply(dat_list, function(df){
    list(seq_len(ncol(df)))
  })
}

# --------------------------------------------------
# Mask time domains by selected locations
# --------------------------------------------------
# domains_by_loc: list(num_locs) of logical(t_length)
# xi_hat: logical(num_locs)
# --------------------------------------------------
# Mask per-location domain selections by a location-selection vector xi
# - domains_by_loc: list length L, each element length T (0/1, TRUE/FALSE, or NA)
# - xi: logical/integer length L (0/1) indicating whether location is selected
# Returns: list length L, each element length T logical (or NA if original was all NA)
mask_domains_by_xi <- function(domains_by_loc, xi, t_length = NULL) {
  stopifnot(is.list(domains_by_loc))
  L <- length(domains_by_loc)
  stopifnot(length(xi) == L)
  
  # infer T if not provided
  if (is.null(t_length)) {
    tt <- vapply(domains_by_loc, length, integer(1))
    # allow some elements to be length 0/NA; choose max nonzero
    t_length <- max(tt, 0L)
  }
  
  out <- vector("list", L)
  
  for (l in seq_len(L)) {
    dom <- domains_by_loc[[l]]
    
    # If domain is NULL, treat as all FALSE (no selections)
    if (is.null(dom)) dom <- rep(FALSE, t_length)
    
    # Coerce to length T
    dom <- as.vector(dom)
    
    # If dom is all NA (your FPCA case), keep as NA mask
    if (all(is.na(dom))) {
      out[[l]] <- rep(NA, t_length)
      next
    }
    
    # If dom length mismatches, pad/truncate safely
    if (length(dom) < t_length) dom <- c(dom, rep(FALSE, t_length - length(dom)))
    if (length(dom) > t_length) dom <- dom[seq_len(t_length)]
    
    dom <- as.logical(dom)
    
    # Apply masking:
    # if xi[l] is FALSE, force all FALSE (no time selected)
    if (!isTRUE(as.logical(xi[l]))) {
      out[[l]] <- rep(FALSE, t_length)
    } else {
      out[[l]] <- dom
    }
  }
  
  out
}

# -------------------------------
# Run one replicate: methods + tests
# -------------------------------

#' Run one replicate and return all key artifacts
#'
#' Methods compared (minimal core set):
#'   - full GP test using all timepoints, all locations
#'   - truncated-by-FocusGP GP test (integrated)
#'   - ISGM selection (on truncated) + optional downstream tests
#'   - ISGM selection (on full) + optional downstream tests
#'
#' @export
run_one_replicate <- function(cfg, demo_stats, real_process = NULL, retry_focusgp = TRUE, max_retries = 5){

  stopifnot(cfg$sim_setting %in% c("real_or_noise", "functional_linear"))
  
  # ----------------------------
  # 0) Known time truth 
  # ----------------------------
  # You can move this into cfg if you want (cfg$true_domain_idx).
  true_domain_idx <- if(!is.null(cfg$true_domain_idx)) cfg$true_domain_idx else 61:80
  
  # 1) simulate dataset
  sim <- simulate_one_dataset(
    cfg = cfg,
    sim_mean_pat = demo_stats$sim_mean_pat,
    sim_sd_pat   = demo_stats$sim_sd_pat,
    sim_mean_ctr = demo_stats$sim_mean_ctr,
    sim_sd_ctr   = demo_stats$sim_sd_ctr,
    real_process = real_process
  )

  full_dat2   <- sim$full_dat2
  sim_scores  <- sim$sim_scores
  true_locs   <- sim$true_locs
  chosen_locs <- sim$chosen_locs
  t_length    <- cfg$t_length
  num_locs    <- cfg$num_locs

  # ----------------------------
  # 1b) Build TRUE DOMAIN(S)
  # ----------------------------
  # Truth domain per location:
  #   - chosen_locs have signal on 61:80 when real_process==1
  #   - other locs have no signal
  true_domains_by_loc <- vector("list", num_locs)
  for(l in seq_len(num_locs)) true_domains_by_loc[[l]] <- rep(FALSE, t_length)
  
  if(cfg$sim_setting == "functional_linear" && sim$real_process == 1){
    for(l in chosen_locs){
      true_domains_by_loc[[l]][true_domain_idx] <- TRUE
    }
  }
  
  # "main" domain truth: evaluate domain selection on FIRST chosen location (matches old selected_time_1 logic)
  if (length(chosen_locs) > 0) {
    domain_loc_eval <- chosen_locs[1]
    true_domains_main <- true_domains_by_loc[[domain_loc_eval]]
  } else {
    domain_loc_eval <- NA_integer_
    true_domains_main <- rep(FALSE, t_length)
  }
  
  # 2) domain selection (FocusGP), with optional retry
  if(retry_focusgp){
    ok <- FALSE
    attempt <- 0
    last_err <- NULL
    while(!ok && attempt < max_retries){
      attempt <- attempt + 1
      out <- tryCatch({
        select_domains_focusgp(full_dat2, sim_scores, cfg)
      }, error = function(e){
        last_err <<- e
        NULL
      })
      if(!is.null(out)) ok <- TRUE
    }
    if(!ok) stop("FocusGP failed after retries: ", last_err$message)
    domsel <- out
  } else {
    domsel <- select_domains_focusgp(full_dat2, sim_scores, cfg)
  }

  # IMPORTANT:
  # domsel$domains is typically a LIST length num_locs (each a length t_length 0/1 vector).
  # For your table replication, we also expose a *single* domain vector for the eval location.
  domains_by_loc <- domsel$domains
  domains_main   <- domains_by_loc[[domain_loc_eval]]
  
  # 3) Build W matrices
  W <- build_W_matrix(
    mode = "simulation",
    full_dat2 = full_dat2,
    chosen_locs = sim$chosen_locs,
    t_length = cfg$t_length
  )
  
  prior_adj <- prior_adjustment(cfg$num_locs, cfg$nx)
  
  # ----------------------------
  # 4) Integrated designs (MATCH DEBUG CODE)
  # ----------------------------
  
  # Truncated: rebuild from full_dat2 using domsel$grps_vec (EXCLUDES intrcpts)
  final_data <- vector("list", num_locs)
  for (l in seq_len(num_locs)) {
    df_l <- as.data.frame(full_dat2[[l]])
    
    # keep exactly the selected time indices from domain selection
    idx <- unname(unlist(domsel$grps_vec[[l]]))
    final_data[[l]] <- df_l[, idx, drop = FALSE]
  }
  
  final_int_df <- as.data.frame(do.call(cbind, final_data))
  grps_vec     <- domsel$grps_vec
  
  # sanity check (this is what was failing)
  map_tr <- build_loc_slices_from_grps(grps_vec)
  stopifnot(sum(map_tr$lens) == ncol(final_int_df))
  
  # Full: numeric-safe wide + matching group format
  final_data_full <- lapply(full_dat2, function(df) {
    as.data.frame(df)[, seq_len(t_length), drop = FALSE]
  })
  final_int_df_full <- as.data.frame(do.call(cbind, final_data_full))
  
  grps_vec_full <- lapply(final_data_full, function(df) {
    list(seq_len(ncol(df)))
  })
  
  map_fl <- build_loc_slices_from_grps(grps_vec_full)
  stopifnot(sum(map_fl$lens) == ncol(final_int_df_full))
  
  # ----------------------------
  # 5) ISGM selection (now consistent)
  # ----------------------------
  curr_xi0 <- rep(1, num_locs)
  
  xi_chain_trunc <- ISGM(
    path_features = final_int_df,
    W_mat = prior_adj * W$W_func,
    y = sim_scores,
    curr_xi = curr_xi0,
    grps_vec = grps_vec,
    connec = cfg$connec,
    N = cfg$isgm_N,
    VBprop = 1
  )
  
  xi_chain_full <- ISGM(
    path_features = final_int_df_full,
    W_mat = prior_adj * W$W_func_full,
    y = sim_scores,
    curr_xi = curr_xi0,
    grps_vec = grps_vec_full,
    connec = cfg$connec,
    N = cfg$isgm_N,
    VBprop = 1
  )

  # 6) Tau tests
  # Full GP design (your original did normalize on as.data.frame(full_dat2))
  X_s <- matrix(1, nrow = cfg$nx, ncol = 1)
  modl_s <- lm(sim_scores ~ X_s - 1)
  GP_df_y <- as.numeric(sim_scores - modl_s$fitted.values)
  
  GP_df_X_full <- apply(as.matrix(final_int_df_full), 2, normalize)
  
  rho_val <- if (is.null(cfg$lengthscale)) {
    estimate_rho_vec(GP_df_X_full)
  } else {
    cfg$lengthscale
  }
  
  p_full <- run_tau_test_fixed(GP_df_X_full, GP_df_y, sigma_vec = estimate_sigma2(GP_df_y), rho = rho_val, X = X_s)

  # Truncated integrated
  GP_df_X_trunc <- as.matrix(apply(final_int_df, 2, normalize))
  
  rho_val <- if (is.null(cfg$lengthscale)) {
    estimate_rho_vec(GP_df_X_trunc)
  } else {
    cfg$lengthscale
  }
  
  p_trunc <- run_tau_test_fixed(GP_df_X_trunc, GP_df_y, sigma_vec = estimate_sigma2(GP_df_y), rho = rho_val, X = X_s)

  # ----------------------------
  # ISGM downstream tau test (TRUNCATED)
  # ----------------------------
  xi_hat_trunc <- as.logical(tail(xi_chain_trunc, 1))  # length = num_locs
  
  # Build column → location map
  block_sizes <- sapply(domsel$grps_vec, function(g) length(unlist(g)))
  block_idx <- rep(seq_along(block_sizes), times = block_sizes)
  
  stopifnot(length(block_idx) == ncol(final_int_df))
  stopifnot(length(xi_hat_trunc) == length(block_sizes))
  
  p_isgm_trunc <- NA_real_
  
  if (sum(xi_hat_trunc) > 0) {
    X_sel <- final_int_df[, block_idx %in% which(xi_hat_trunc), drop = FALSE]
    X_sel <- apply(as.matrix(X_sel), 2, normalize)
    
    rho_val <- if (is.null(cfg$lengthscale)) {
      estimate_rho_vec(X_sel)
    } else {
      cfg$lengthscale
    }
    
    p_isgm_trunc <- run_tau_test_fixed(
      X_sel,
      GP_df_y,
      sigma_vec = estimate_sigma2(GP_df_y),
      rho = rho_val,
      X = X_s
    )
  }
  
  xi_hat_full <- as.logical(tail(xi_chain_full, 1))
  
  block_sizes_full <- rep(cfg$t_length, cfg$num_locs)
  block_idx_full <- rep(seq_len(cfg$num_locs), times = block_sizes_full)
  
  p_isgm_full <- NA_real_
  
  if (sum(xi_hat_full) > 0) {
    X_sel <- final_int_df_full[, block_idx_full %in% which(xi_hat_full), drop = FALSE]
    X_sel <- apply(as.matrix(X_sel), 2, normalize)
    
    rho_val <- if (is.null(cfg$lengthscale)) {
      estimate_rho_vec(X_sel)
    } else {
      cfg$lengthscale
    }
    
    p_isgm_full <- run_tau_test_fixed(
      X_sel,
      GP_df_y,
      sigma_vec = estimate_sigma2(GP_df_y),
      rho = rho_val,
      X = X_s
    )
  }
  
  # --- your method has domains per location (from FocusGP) ---
  # Create “flavors”:
  domains_by_loc_raw        <- domains_by_loc
  domains_by_loc_isgm_trunc <- mask_domains_by_xi(domains_by_loc, xi_hat_trunc, t_length = t_length)
  domains_by_loc_isgm_full <- vector("list", num_locs)
  for (l in seq_len(num_locs)) {
    if (isTRUE(xi_hat_full[l])) {
      domains_by_loc_isgm_full[[l]] <- rep(TRUE, t_length)
    } else {
      domains_by_loc_isgm_full[[l]] <- rep(FALSE, t_length)
    }
  }

  # 7) Return a single object per replicate
  list(
    meta = list(
      real_process = sim$real_process,
      true_locs = true_locs,
      chosen_locs = sim$chosen_locs,
      prior_adj = prior_adj
    ),
    full_dat2   = full_dat2,
    sim_scores  = sim_scores,
    
    # keep original if you want
    domains_by_loc = domains_by_loc,
    
    # NEW: masked variants (these are what you will pass to log_domain_method_per_loc)
    domains_by_loc_raw        = domains_by_loc_raw,
    domains_by_loc_isgm_trunc = domains_by_loc_isgm_trunc,
    domains_by_loc_isgm_full  = domains_by_loc_isgm_full,
    
    xi_chain_trunc = xi_chain_trunc,
    xi_chain_full  = xi_chain_full,
    xi_hat_trunc   = xi_hat_trunc,
    xi_hat_full    = xi_hat_full,
    tau_pvals = list(
      full       = p_full,
      trunc      = p_trunc,
      isgm_trunc = p_isgm_trunc,
      isgm_full  = p_isgm_full
    )
  )
}


# -------------------------------
# Monte Carlo runner
# -------------------------------

#' Run a Monte Carlo simulation
#'
#' @param cfg config
#' @param demo_stats list with sim_mean_pat, sim_sd_pat, sim_mean_ctr, sim_sd_ctr
#' @param realp_vec optional vector length sim_iters (0/1); if NULL, random
#' @export
run_mc <- function(cfg, demo_stats, realp_vec = NULL, verbose = TRUE){

  S <- cfg$sim_iters
  out <- vector("list", length = S)

  for(s in seq_len(S)){
    if(verbose) message("Simulation iteration ", s, "/", S)

    rp <- NULL
    if(!is.null(realp_vec)) rp <- realp_vec[s]

    out[[s]] <- run_one_replicate(
      cfg = cfg,
      demo_stats = demo_stats,
      real_process = rp,
      retry_focusgp = TRUE,
      max_retries = 5
    )
  }

  out
}


# -------------------------------
# Method comparison + simple stats
# -------------------------------

#' Summarize p-values and selection accuracy across MC reps
#'
#' Produces:
#'  - per-method power (if real_process==1)
#'  - per-method type I error (if real_process==0)
#'  - location selection accuracy (TP/FP counts)
#'
#' @export
summarize_mc <- function(mc_out, alpha = 0.05){

  # helper to coerce tau p-value output
  get_p <- function(p){
    # your tau.score.test currently returns a p-value (scalar). Some variants return list.
    if(is.list(p) && !is.null(p$p.value)) return(as.numeric(p$p.value))
    as.numeric(p)
  }

  S <- length(mc_out)

  methods <- c("full", "trunc", "isgm_trunc", "isgm_full")
  pmat <- matrix(NA_real_, nrow = S, ncol = length(methods))
  colnames(pmat) <- methods

  realp <- integer(S)
  tp_trunc <- fp_trunc <- tp_full <- fp_full <- integer(S)

  for(i in seq_len(S)){
    realp[i] <- mc_out[[i]]$meta$real_process

    for(m in methods){
      pmat[i, m] <- get_p(mc_out[[i]]$tau_pvals[[m]])
    }

    true <- mc_out[[i]]$meta$true_locs

    x_tr <- mc_out[[i]]$xi_hat_trunc
    x_fl <- mc_out[[i]]$xi_hat_full

    tp_trunc[i] <- sum((x_tr == 1) & (true == 1))
    fp_trunc[i] <- sum((x_tr == 1) & (true == 0))
    tp_full[i]  <- sum((x_fl == 1) & (true == 1))
    fp_full[i]  <- sum((x_fl == 1) & (true == 0))
  }

  # power / type I error
  power <- sapply(methods, function(m){
    mean(pmat[realp == 1, m] < alpha, na.rm = TRUE)
  })

  type1 <- sapply(methods, function(m){
    mean(pmat[realp == 0, m] < alpha, na.rm = TRUE)
  })

  list(
    meta = list(
      real_process = sim$real_process,
      true_locs    = true_locs,
      chosen_locs  = sim$chosen_locs,
      prior_adj    = prior_adj
    ),
    full_dat2  = full_dat2,
    sim_scores = sim_scores,
    
    # Domain outputs
    domains_by_loc        = domains_by_loc,          # raw FocusGP
    domains_by_loc_masked = domains_by_loc_masked,   # ✅ FIXED output
    domains_main          = domains_main,             # for legacy comparison
    
    # Node outputs
    xi_hat_trunc = xi_hat_trunc,
    xi_hat_full  = xi_hat_full,
    
    tau_pvals = list(
      full       = p_full,
      trunc      = p_trunc,
      isgm_trunc = p_isgm_trunc,
      isgm_full  = p_isgm_full
    )
  )
}
