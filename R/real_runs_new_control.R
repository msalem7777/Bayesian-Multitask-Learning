rm(list = ls())
graphics.off()

source("C:/My Files/repos/Bayesian-Multitask-Learning/R/load_libraries.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/tau_tests.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/run_simulation.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/data_generation.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/domain_selection.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/domain_selection_core.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/correlation_matrix_constructor.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/location_selection.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/ISGM_core.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/evaluation_metrics.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/baseline_functions.R")

# --------------------------------------------------
# Config
# --------------------------------------------------

brain_regions <- read.csv("C:/My files/Transfer files/PhD Research/P2/Data/Set of FMRI/Network/Documents 2/Schaefer2018_100Parcels_7Networks_2mm.csv", header = TRUE)
brain_regions_unique = unique(brain_regions$Network.Name)

score_column    <- "ADOS_Total"
t_length        <- 100
sim_iters       <- 1    # bootstrap-repeated analysis runs; set to 1 for a single full-cohort run
isgm_N          <- 300   # ISGM MCMC iterations
connec          <- 1.0   # ISGM connectivity scaling


# --------------------------------------------------
# Load ABIDE demographics and filter
# --------------------------------------------------

dat_ABIDE_demo <- read.csv(
  "C:/My Files/repos/Bayesian-Multitask-Learning/R/Demographics_ABIDE1_2_extended_F4.csv"
)

dat_clipped <- dat_ABIDE_demo[, c(
  "subID", "SubjectType", "AgeAtScan", "TRSec",
  "SamplingRate", "Sex", "Handedness_Category",
  "ADOS_Total", "SRS_RawTotal"
)]

ADOS_dat   <- dat_clipped[!is.na(dat_clipped$ADOS_Total), ]
dat_active <- ADOS_dat[ADOS_dat$SamplingRate == 0.5, ]

dat_active$SubjectType <- as.factor(dat_active$SubjectType)
dat_active$Sex         <- as.factor(dat_active$Sex)
dat_active$Handedness  <- as.factor(dat_active$Handedness_Category)

N     <- length(unique(dat_active$subID))
N_pat <- sum(dat_active$SubjectType == "PATIENT")
N_ctr <- sum(dat_active$SubjectType == "CONTROL")

for(br_reg in brain_regions_unique){

  brain_selection <- br_reg
  if (brain_selection == "Salience / Ventral Attention"){
    brain_write = "SV"
  } else{
    brain_write = brain_selection
  }
  
  # --------------------------------------------------
  # Load brain-region list and fMRI data
  # --------------------------------------------------
  
  brain_regions_selected <- brain_regions$ROI.Label[
    brain_regions$Network.Name == brain_selection
  ]
  num_locs <- length(brain_regions_selected)
  
  # dat_listr_AAL2_tilde[[l]]: data frame, rows = subjects, columns:
  #   ID, t1..t100, status, sex, age, handedness, score
  #   time columns have been residualised for sex / age / handedness
  dat_listr_AAL2       <- list()
  dat_listr_AAL2_tilde <- list()
  dat_listr_coeffs     <- list()
  
  for (l in brain_regions_selected) {
    
    dat_holdr <- data.frame(matrix(NA, nrow = N, ncol = t_length + 2))
    colnames(dat_holdr) <- c("ID", paste0("t", seq_len(t_length)), "status")
    dat_holdr$ID         <- dat_active$subID
    dat_holdr$status     <- dat_active$SubjectType
    dat_holdr$sex        <- dat_active$Sex
    dat_holdr$age        <- dat_active$AgeAtScan
    dat_holdr$handedness <- dat_active$Handedness
    dat_holdr$score      <- dat_active[[score_column]]
    
    for (i in seq_len(N)) {
      fmri_path <- paste0(
        "C:/My files/Transfer files/PhD Research/P2/Data/Set of FMRI/FMRI/Data_ABIDE1_2_extended_csv/Ind_data_extended/",
        dat_active$subID[i], "/Schaefer100.csv"
      )
      
      fMRI_dat <- read.csv(fmri_path, header = FALSE)
      dat_holdr[i, 2:(t_length + 1)] <- fMRI_dat[1:t_length, l]
    }
    
    # Regress out sex, age, handedness (and sex × age interaction)
    dat_holdr_long <- melt(
      setDT(dat_holdr),
      id.vars       = c("ID", "status", "sex", "age", "handedness", "score"),
      variable.name = "time"
    )
    model <- lm(value ~ sex + age + handedness + sex * age, data = dat_holdr_long)
    cofs  <- coefficients(model)
    
    dat_holdr_tilde <- dat_holdr
    for (i in seq_len(N)) {
      dat_holdr_tilde[i, 2:(t_length + 1)] <-
        dat_holdr[i, 2:(t_length + 1)] -
        rep((cofs["sex2"])            * (dat_holdr[i, "sex"] == 2),                         t_length) -
        rep((cofs["sex2:age"])        * (dat_holdr[i, "sex"] == 2) * dat_holdr[i, "age"],   t_length) -
        rep((cofs["handednessR"])     * (dat_holdr[i, "handedness"] == "R"),                t_length) -
        rep((cofs["handednessMixed"]) * (dat_holdr[i, "handedness"] == "Mixed"),            t_length) -
        rep((cofs["age"])             * dat_holdr[i, "age"],                                t_length)
    }
    
    dat_listr_coeffs     <- c(dat_listr_coeffs, list(cofs))
    dat_listr_AAL2       <- c(dat_listr_AAL2,   list(dat_holdr))
    dat_listr_AAL2_tilde <- c(dat_listr_AAL2_tilde, list(dat_holdr_tilde))
  }
  
  
  # --------------------------------------------------
  # Helper: extract full_dat2 from the tilde list.
  #
  # run_simulation.R convention: full_dat2[[l]] is [N × t_length] (subjects as rows),
  # matching dat_listr_AAL2_tilde[[l]][, 2:(t_length+1)].
  # No transposition needed.
  # --------------------------------------------------
  
  make_full_dat2 <- function(dat_listr_tilde, num_locs, t_length) {
    lapply(seq_len(num_locs), function(l) {
      as.data.frame(dat_listr_tilde[[l]][, 2:(t_length + 1)])   # [N × t_length]
    })
  }
  
  
  # --------------------------------------------------
  # Helper: residualise outcome and build projection matrix.
  # Mirrors run_one_replicate()'s tau-test setup, but uses
  # real covariates (age, sex, handedness) in place of
  # the simulation's trivial intercept-only X_s.
  # --------------------------------------------------
  
  residualise_scores <- function(dat_active, sim_scores) {
    
    df_dat <- data.frame(
      score       = sim_scores,
      age         = dat_active$AgeAtScan,
      sex         = dat_active$Sex,
      handedness  = dat_active$Handedness
    )
    modl_s  <- lm(score ~ age + sex + handedness, data = df_dat)
    y_resid <- as.numeric(sim_scores - modl_s$fitted.values)
    
    # Design matrix for tau.score.test's X argument (projection-out of confounders)
    X_s <- model.matrix(~ age + sex + handedness, data = df_dat)
    
    sigma_vec <- as.numeric((1 / length(y_resid)) * crossprod(y_resid))
    
    list(y_resid = y_resid, sigma_vec = sigma_vec, X_s = X_s)
  }
  
  
  # --------------------------------------------------
  # Core per-replicate function — mirrors run_one_replicate()
  # --------------------------------------------------
  
  run_one_real_replicate <- function(full_dat2, dat_active, sim_scores,
                                     num_locs, t_length, N, W_mat, isgm_N, connec) {
    
    # ---- 2. Domain selection via select_domains_focusgp() ----
    # Reuse the same function as run_one_replicate(); it reads full_dat2[[l]][, 1:t_length].
    # We wrap it in the same retry logic.
    cfg_domain <- list(
      num_locs       = num_locs,
      t_length       = t_length,
      nx             = N,
      focusgp_nfold  = 5,
      focusgp_lambda = 0.05,
      focusgp_ninit  = 12,
      focusgp_end    = 30
    )
    
    domsel <- NULL
    for (attempt in seq_len(5)) {
      domsel <- tryCatch(
        select_domains_focusgp(full_dat2, sim_scores, cfg_domain),
        error = function(e) NULL
      )
      if (!is.null(domsel)) break
    }
    if (is.null(domsel)) stop("FocusGP failed after 5 retries")
    
    domains_by_loc_raw <- domsel$domains   # list(num_locs), each length-t_length 0/1
    grps_vec           <- domsel$grps_vec  # list(num_locs) of consecutive-index groups
    
    # ---- 3. Build integrated design matrices (mirrors run_one_replicate step 4) ----
    # Truncated: time indices selected by FocusGP per location
    final_data_trunc <- vector("list", num_locs)
    for (l in seq_len(num_locs)) {
      idx <- unname(unlist(grps_vec[[l]]))
      final_data_trunc[[l]] <- as.data.frame(full_dat2[[l]])[, idx, drop = FALSE]
    }
    final_int_df <- as.data.frame(do.call(cbind, final_data_trunc))
    
    # Full: all timepoints
    final_data_full <- lapply(full_dat2, function(df) {
      as.data.frame(df)[, seq_len(t_length), drop = FALSE]
    })
    final_int_df_full <- as.data.frame(do.call(cbind, final_data_full))
    
    # grps_vec variants for ISGM
    grps_vec_full <- lapply(final_data_full, function(df) list(seq_len(ncol(df))))
    
    # ---- 4. ISGM location selection (truncated + full), matching run_one_replicate step 5 ----
    prior_adj <- prior_adjustment(num_locs, N)
    curr_xi0  <- rep(1L, num_locs)
    
    xi_chain_trunc <- ISGM(
      path_features = final_int_df,
      W_mat         = prior_adj * W_mat,
      y             = sim_scores,
      curr_xi       = curr_xi0,
      grps_vec      = grps_vec,
      connec        = connec,
      N             = isgm_N,
      VBprop        = 1
    )
    
    xi_chain_full <- ISGM(
      path_features = final_int_df_full,
      W_mat         = prior_adj * W_mat,
      y             = sim_scores,
      curr_xi       = curr_xi0,
      grps_vec      = grps_vec_full,
      connec        = connec,
      N             = isgm_N,
      VBprop        = 1
    )
    
    xi_hat_trunc <- as.logical(tail(xi_chain_trunc, 1))
    xi_hat_full  <- as.logical(tail(xi_chain_full,  1))
    
    # ---- 5. Domain mask variants (mirrors run_one_replicate step 6 end) ----
    domains_by_loc_isgm_trunc <- mask_domains_by_xi(domains_by_loc_raw, xi_hat_trunc, t_length)
    
    domains_by_loc_isgm_full <- vector("list", num_locs)
    for (l in seq_len(num_locs)) {
      domains_by_loc_isgm_full[[l]] <- if (isTRUE(xi_hat_full[l])) rep(TRUE, t_length) else rep(FALSE, t_length)
    }
    
    # ---- 6. Tau score tests with covariate residualisation ----
    resid_out <- residualise_scores(dat_active, sim_scores)
    y_resid   <- resid_out$y_resid
    sigma_vec <- resid_out$sigma_vec
    X_s       <- resid_out$X_s
    
    # Mirrors run_one_replicate's block-index + estimate_rho_vec pattern
    run_tau <- function(feature_df, block_idx, xi_sel) {
      if (sum(xi_sel) == 0) return(NA_real_)
      X_sel <- feature_df[, block_idx %in% which(xi_sel), drop = FALSE]
      X_sel <- apply(as.matrix(X_sel), 2, normalize)
      rho   <- 4 #estimate_rho_vec(X_sel)
      print(rho)
      run_tau_test_fixed(X_sel, y_resid, sigma_vec = sigma_vec, rho = rho, X = X_s)
    }
    
    # Build block-index vectors
    block_sizes_trunc <- sapply(grps_vec, function(g) length(unlist(g)))
    block_idx_trunc   <- rep(seq_along(block_sizes_trunc), times = block_sizes_trunc)
    
    block_idx_full <- rep(seq_len(num_locs), each = t_length)
    
    tau_full       <- run_tau(final_int_df_full, block_idx_full,  rep(TRUE, num_locs))
    tau_trunc      <- run_tau(final_int_df,      block_idx_trunc, rep(TRUE, num_locs))
    tau_isgm_trunc <- run_tau(final_int_df,      block_idx_trunc, xi_hat_trunc)
    tau_isgm_full  <- run_tau(final_int_df_full, block_idx_full,  xi_hat_full)
    
    # ---- Return: structure mirrors run_one_replicate() ----
    list(
      full_dat2                 = full_dat2,
      sim_scores                = sim_scores,
      domains_by_loc_raw        = domains_by_loc_raw,
      domains_by_loc_isgm_trunc = domains_by_loc_isgm_trunc,
      domains_by_loc_isgm_full  = domains_by_loc_isgm_full,
      xi_hat_trunc              = xi_hat_trunc,
      xi_hat_full               = xi_hat_full,
      tau_pvals = list(
        full       = tau_full,
        trunc      = tau_trunc,
        isgm_trunc = tau_isgm_trunc,
        isgm_full  = tau_isgm_full
      )
    )
  }
  
  
  # --------------------------------------------------
  # Replicate loop
  # --------------------------------------------------
  
  sim_scores <- dat_active[[score_column]]
  
  domain_rows <- list()
  node_rows   <- list()
  k_domain    <- 1
  k_node      <- 1
  
  # Null placeholders — no ground truth for real data
  null_truth_domains <- replicate(num_locs, rep(NA, t_length), simplify = FALSE)
  null_true_locs     <- rep(NA, num_locs)
  null_meta          <- list(
    num_locs    = num_locs,
    t_length    = t_length,
    chosen_locs = integer(0),
    real_process = NA
  )
  
  for (s in seq_len(sim_iters)) {
    
    # set.seed(s)
    # cat("--- Real data replicate", s, "of", sim_iters, "---\n")
    
    # Bootstrap resample for repeated-analysis stability assessment.
    # Remove the next three lines and set dat_active_s/dat_listr_s/sim_scores_s
    # directly from the full-cohort objects for a single non-resampled run.
    # bs_idx       <- sample(seq_len(N), N, replace = TRUE)
    # dat_active_s <- dat_active[bs_idx, ]
    # dat_listr_s  <- lapply(dat_listr_AAL2_tilde, function(df) df[bs_idx, ])
    # sim_scores_s <- sim_scores[bs_idx]
    
    dat_active_s <- dat_active
    dat_listr_s  <- dat_listr_AAL2_tilde
    sim_scores_s <- sim_scores
    
    full_dat2_s <- make_full_dat2(dat_listr_s, num_locs, t_length)
    
    W_mat <- build_W_matrix(
      mode      = "data",
      full_dat2 = full_dat2_s,
      t_length  = t_length
    )$W_func
    
    res <- run_one_real_replicate(
      full_dat2  = full_dat2_s,
      dat_active = dat_active_s,
      sim_scores = sim_scores_s,
      num_locs   = num_locs,
      t_length   = t_length,
      N          = N,
      W_mat      = W_mat,
      isgm_N     = isgm_N,
      connec     = connec
    )
    
    rep_id <- paste0("R", s)
    
    # ---- Domain rows ----
    for (meth_pair in list(
      list("DomainSelection_FocusGP_RAW",        res$domains_by_loc_raw),
      list("DomainSelection_FocusGP_ISGM_TRUNC", res$domains_by_loc_isgm_trunc),
      list("DomainSelection_FocusGP_ISGM_FULL",  res$domains_by_loc_isgm_full)
    )) {
      domain_rows[[k_domain]] <- log_domain_method_per_loc(
        method                  = meth_pair[[1]],
        selected_domains_by_loc = meth_pair[[2]],
        truth_domains_by_loc    = null_truth_domains,
        meta                    = null_meta,
        rep_id                  = rep_id
      )
      k_domain <- k_domain + 1
    }
    
    # Baselines
    lasso_out <- run_lasso(res$full_dat2, res$sim_scores, num_locs, t_length)
    fl_out    <- run_fused_lasso(res$full_dat2, res$sim_scores, num_locs)
    bs13_out  <- run_bspline_lm(res$full_dat2, res$sim_scores, num_locs, t_length, bn = 13)
    bs23_out  <- run_bspline_lm(res$full_dat2, res$sim_scores, num_locs, t_length, bn = 23)
    fpca_locs <- run_gp_fpca(res$full_dat2, res$sim_scores, num_locs, t_length, N)
    
    for (meth_pair in list(
      list("Lasso",      lasso_out$domains_by_loc),
      list("FusedLasso", fl_out$domains_by_loc),
      list("FLM_BS13",   bs13_out$domains_by_loc),
      list("FLM_BS23",   bs23_out$domains_by_loc),
      list("FPCA_GP",    replicate(num_locs, rep(NA, t_length), simplify = FALSE))
    )) {
      domain_rows[[k_domain]] <- log_domain_method_per_loc(
        method                  = meth_pair[[1]],
        selected_domains_by_loc = meth_pair[[2]],
        truth_domains_by_loc    = null_truth_domains,
        meta                    = null_meta,
        rep_id                  = rep_id
      )
      k_domain <- k_domain + 1
    }
    
    # ---- Node rows ----
    node_rows[[k_node]] <- log_method(
      method        = "A1_Full_GP",
      selected_locs = rep(TRUE, num_locs),
      tau           = res$tau_pvals$full,
      truth         = null_true_locs,
      meta          = null_meta,
      rep_id        = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      method        = "A2_DomainOnly",
      selected_locs = rep(TRUE, num_locs),
      tau           = res$tau_pvals$trunc,
      truth         = null_true_locs,
      meta          = null_meta,
      rep_id        = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      method        = "A3_Domain_ISGM",
      selected_locs = res$xi_hat_trunc,
      tau           = res$tau_pvals$isgm_trunc,
      truth         = null_true_locs,
      meta          = null_meta,
      rep_id        = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      method        = "A4_ISGM_Full",
      selected_locs = res$xi_hat_full,
      tau           = res$tau_pvals$isgm_full,
      truth         = null_true_locs,
      meta          = null_meta,
      rep_id        = rep_id
    )
    k_node <- k_node + 1
    
    for (meth_pair in list(
      list("Lasso",      lasso_out$selected_locs),
      list("FusedLasso", fl_out$selected_locs),
      list("BSpline_13", bs13_out$selected_locs),
      list("BSpline_23", bs23_out$selected_locs),
      list("GP_FPCA",    fpca_locs)
    )) {
      node_rows[[k_node]] <- log_method(
        method        = meth_pair[[1]],
        selected_locs = meth_pair[[2]],
        tau           = NA_real_,
        truth         = null_true_locs,
        meta          = null_meta,
        rep_id        = rep_id
      )
      k_node <- k_node + 1
    }
  }
  
  
  # --------------------------------------------------
  # Collect and save
  # --------------------------------------------------
  
  domain_results_df <- do.call(rbind, domain_rows)
  node_results_df   <- do.call(rbind, node_rows)
  
  print(domain_results_df)
  print(node_results_df)
  
  
  out_path <- file.path(
    "C:/My Files/repos/Bayesian-Multitask-Learning/R",
    paste0("real_domain_results_", brain_write, ".csv")
  )
  
  write.csv(domain_results_df, out_path, row.names = FALSE)
  
  
  out_path <- file.path(
    "C:/My Files/repos/Bayesian-Multitask-Learning/R",
    paste0("real_node_results_", brain_write, ".csv")
  )
  
  write.csv(node_results_df, out_path, row.names = FALSE)
  
  
  # ---- Export summary table (6 rows × num_locs columns) ----
  
  domain_to_range <- function(dom) {
    idx <- which(as.logical(dom))
    if (length(idx) == 0) return("")
    paste0("(", min(idx), ",", max(idx), ")")
  }
  
  loc_names <- brain_regions$ROI.Name[brain_regions_selected]
  
  xi_domain_export <- data.frame(
    rbind(
      as.character(as.integer(res$xi_hat_trunc)),
      as.character(as.integer(res$xi_hat_full)),
      sapply(res$domains_by_loc_isgm_trunc, domain_to_range),
      rep(as.character(res$tau_pvals$isgm_trunc), num_locs),
      rep(as.character(res$tau_pvals$isgm_full),  num_locs)
    ),
    stringsAsFactors = FALSE
  )
  colnames(xi_domain_export)  <- loc_names
  rownames(xi_domain_export)  <- c("selected_trunc", "selected_full", "domain_trunc", "pval_trunc", "pval_full")
  
  out_path <- file.path(
    "C:/My Files/repos/Bayesian-Multitask-Learning/R",
    paste0("real_xi_domain_export_", brain_write, ".csv")
  )
  
  write.csv(xi_domain_export, out_path, row.names = TRUE)
  
  # ---- GP residual plot ----
  resid_out <- residualise_scores(dat_active_s, sim_scores_s)
  y_resid   <- resid_out$y_resid
  plot(y_resid, col = dat_active_s$SubjectType)
  
  X_sel <- do.call(cbind, lapply(which(res$xi_hat_trunc), function(l) {
    idx <- which(as.logical(res$domains_by_loc_isgm_trunc[[l]]))
    as.matrix(res$full_dat2[[l]])[, idx, drop = FALSE]
  }))
  X_sel <- apply(normalize(X_sel), 2, normalize)
  
  gpi <- newGPsep(X_sel, y_resid, d = 4, g = 0.001 * var(y_resid), dK = TRUE)
  # mleGPsep(gpi, param = "d", tmin = 0.1, tmax = 1e4)
  gp_pred  <- predGPsep(gpi, X_sel)$mean
  deleteGPsep(gpi)
  
  gp_resid <- y_resid - gp_pred
  
  plot(
    gp_pred, gp_resid,
    xlab = "GP fitted values",
    ylab = "GP residuals",
    main = paste0("GP residual plot — ", brain_selection, " (replicate ", s, ")"),
    pch  = 19,
    col  = ifelse(dat_active_s$SubjectType == "PATIENT", "tomato", "steelblue")
  )
  abline(h = 0, lty = 2)
  legend("topright", legend = c("Patient", "Control"),
         col = c("tomato", "steelblue"), pch = 19, bty = "n")
}