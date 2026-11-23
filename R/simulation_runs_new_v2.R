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


# Loading the ABIDE AUTISM demographic data
dat_ABIDE_demo = read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/Demographics_ABIDE1_2_extended_F4.csv")
setwd("C:/My Files/repos/Bayesian-Multitask-Learning/")

dat_clipped = dat_ABIDE_demo[,c("subID", "SubjectType", "AgeAtScan", "TRSec","SamplingRate", "Sex", "Handedness_Category", "ADOS_Total", "SRS_RawTotal")]
SRS_dat = dat_clipped[!is.na(dat_clipped$SRS_RawTotal),]
ADOS_dat = dat_clipped[!is.na(dat_clipped$ADOS_Total),]

dat_active = ADOS_dat
dat_active = dat_active[dat_active$SamplingRate==0.5,]

sim_mean_pat = mean(dat_active$ADOS_Total[dat_active$SubjectType=="PATIENT"])
sim_sd_pat = sd(dat_active$ADOS_Total[dat_active$SubjectType=="PATIENT"])
sim_mean_ctr = mean(dat_active$ADOS_Total[dat_active$SubjectType=="CONTROL"])
sim_sd_ctr = sd(dat_active$ADOS_Total[dat_active$SubjectType=="CONTROL"])
# --------------------------------------------------
# Config
# --------------------------------------------------
cfg <- make_default_cfg(
  seed = 1,
  sim_iters = 100,
  sim_setting = "real_or_noise",
  nx = 200,
  num_locs = 10,
  num_T = 20,
  chosen_locs = c(2,5,8),
  lengthscale = 1, 
  chosen_one = 1,
  unchosen_one = 2
)

demo_stats <- list(
  sim_mean_pat = sim_mean_pat,
  sim_sd_pat   = sim_sd_pat,
  sim_mean_ctr = sim_mean_ctr,
  sim_sd_ctr   = sim_sd_ctr
)

t_length = 100
# -------------------------------------------------
# GP-on-FPCA
# --------------------------------------------------
run_gp_fpca <- function(full_dat2, sim_scores, num_locs, t_length, N, thresh = 8.0) {
  
  library(fdapace)
  
  pc_mat <- matrix(0, nrow = N, ncol = 0)
  t_list <- replicate(N, seq_len(t_length), simplify = FALSE)
  
  pc_length <- integer(0)
  
  for (l in seq_len(num_locs)) {
    
    y_list = as.data.frame(t(full_dat2[[l]]))
    y_list = as.list(y_list)
    
    FPCAobj <- FPCA(
      Ly = y_list,
      Lt = t_list,
      optns = list(dataType = "Dense")
    )
    
    pc_mat <- cbind(pc_mat, FPCAobj$xiEst)
    pc_length <- c(pc_length, ncol(FPCAobj$xiEst))
  }

  pc_length_cum = c(0,cumsum(pc_length))

  gpi <- newGPsep(
    normalize(pc_mat),
    as.numeric(scale(sim_scores)),
    d = 0.1,
    g = 0.001,
    dK = TRUE
  )
  
  mle <- mleGPsep(
    gpi,
    param = "both",
    tmin = c(0.0001, 0.000001),
    tmax = c(10, 10)
  )
  selected_locs <- logical(num_locs)
  for (l in seq_len(num_locs)) {
    selected_locs[l] = 1*(min(mle$theta[(pc_length_cum[l]+1):pc_length_cum[l+1]], na.rm=TRUE)<8.0)
  }

  deleteGPseps()
  
  selected_locs
}


# -------------------------------------------------
# Bspline
# --------------------------------------------------
# legacy-like mapping
bspline_basis_to_time_mask <- function(pvals_lb, t_length, bn) {
  mask <- rep(FALSE, t_length)
  
  if (bn == 23) {
    # exactly like your old bn=23 mapping (step=5, window=20, edges special)
    for (lb in seq_len(bn)) {
      if (!(pvals_lb[lb] < 0.05)) next
      
      if (lb == 1) {
        mask[1:5] <- TRUE
      } else if (lb == 2) {
        mask[1:10] <- TRUE
      } else if (lb == 3) {
        mask[1:15] <- TRUE
      } else if (lb == 23) {
        mask[95:100] <- TRUE
      } else if (lb == 22) {
        mask[90:100] <- TRUE
      } else if (lb == 21) {
        mask[85:100] <- TRUE
      } else {
        a <- (lb - 4) * 5 + 1
        b <- (lb - 4) * 5 + 20
        mask[a:b] <- TRUE
      }
    }
    return(mask)
  }
  
  if (bn == 13) {
    # exactly like your old bn=13 mapping (step=10, window=40, edges special)
    for (lb in seq_len(bn)) {
      if (!(pvals_lb[lb] < 0.05)) next
      
      if (lb == 1) {
        mask[1:10] <- TRUE
      } else if (lb == 2) {
        mask[1:20] <- TRUE
      } else if (lb == 3) {
        mask[1:30] <- TRUE
      } else if (lb == 13) {
        mask[90:100] <- TRUE
      } else if (lb == 12) {
        mask[80:100] <- TRUE
      } else if (lb == 11) {
        mask[70:100] <- TRUE
      } else {
        a <- (lb - 4) * 10 + 1
        b <- (lb - 4) * 10 + 40
        mask[a:b] <- TRUE
      }
    }
    return(mask)
  }
  
  stop("Unsupported bn for legacy mapping: ", bn)
}

run_bspline_lm <- function(full_dat2, sim_scores, num_locs, t_length, bn) {
  library(fda)
  
  N <- length(sim_scores)
  
  dp_mat <- matrix(0, nrow = N)
  basisobj <- create.bspline.basis(c(1, t_length), bn)
  
  for (l in seq_len(num_locs)) {
    fd_dat <- t(full_dat2[[l]])
    colnames(fd_dat) <- seq_len(N)
    rownames(fd_dat) <- seq_len(t_length)
    
    b_fd <- smooth.basisPar(seq_len(t_length), fd_dat, basisobj)$fd
    dp_mat <- cbind(dp_mat, t(as.matrix(b_fd$coefs)))
  }
  
  dp_mat <- dp_mat[, 2:ncol(dp_mat)]
  colnames(dp_mat) <- paste0("X", seq_len(ncol(dp_mat)))
  
  fit <- lm(scale(sim_scores) ~ dp_mat)
  fit_summ <- summary(fit)
  
  pvals <- fit_summ$coefficients[2:nrow(fit_summ$coefficients), 4]
  
  # NA handling same as your code
  pvals_nan <- seq_along(pvals) * is.na(pvals)
  pvals[pvals_nan] <- 1
  if (length(pvals) != ncol(dp_mat)) {
    pvals <- c(pvals, rep(1, ncol(dp_mat) - length(pvals)))
  }
  
  selected_locs  <- logical(num_locs)
  domains_by_loc <- vector("list", num_locs)
  
  for (lnum in seq_len(num_locs)) {
    idx <- ((lnum - 1) * bn + 1):(lnum * bn)
    pvals_lb <- pvals[idx]
    
    selected_locs[lnum] <- (min(pvals_lb) < 0.05)
    domains_by_loc[[lnum]] <- bspline_basis_to_time_mask(pvals_lb, t_length, bn)
  }
  
  list(
    selected_locs = selected_locs,
    domains_by_loc = domains_by_loc
  )
}

# --------------------------------------------------
# Fused Lasso
# --------------------------------------------------
run_fused_lasso <- function(full_dat2, sim_scores, num_locs,
                            lambda1 = 0.1, lambda2 = 1, thresh = 0.005) {
  
  library(extlasso)
  
  selected_locs  <- logical(num_locs)
  domains_by_loc <- vector("list", num_locs)
  
  for (l in seq_len(num_locs)) {
    fit <- extlasso::fusedlasso(as.matrix(full_dat2[[l]]),
                                scale(sim_scores),
                                lambda1 = lambda1,
                                lambda2 = lambda2,
                                maxiter = 2e5)
    
    dom <- abs(fit$coef) > thresh
    domains_by_loc[[l]] <- as.logical(dom)
    selected_locs[l] <- any(dom)
  }
  
  list(
    selected_locs = selected_locs,
    domains_by_loc = domains_by_loc
  )
}


# --------------------------------------------------
# Lasso
# --------------------------------------------------
run_lasso <- function(full_dat2, sim_scores, num_locs, t_length) {
  
  library(glmnet)
  
  big_DF <- as.data.frame(do.call(cbind, lapply(full_dat2, as.matrix)))
  
  cv_model <- cv.glmnet(as.matrix(big_DF), scale(sim_scores), alpha = 1)
  best_model <- glmnet(as.matrix(big_DF), scale(sim_scores),
                       alpha = 1, lambda = cv_model$lambda.min)
  
  beta <- rep(0, ncol(big_DF))
  beta[best_model$beta@i] <- best_model$beta@x
  
  # per-location domain masks
  domains_by_loc <- vector("list", num_locs)
  selected_locs  <- logical(num_locs)
  
  for (l in seq_len(num_locs)) {
    idx <- ((l - 1) * t_length + 1):(l * t_length)
    dom <- beta[idx] != 0
    domains_by_loc[[l]] <- as.logical(dom)
    selected_locs[l] <- any(dom)
  }
  
  list(
    selected_locs = selected_locs,
    domains_by_loc = domains_by_loc,
    beta = beta
  )
}

# --------------------------------------------------
# Loop over replicates
# --------------------------------------------------
domain_rows <- list()
node_rows   <- list()

k_domain <- 1
k_node   <- 1

for (scenario in c(0, 1)) {   # 0 = Null, 1 = Alternative

  for (s in seq_len(cfg$sim_iters)) {
    
    cfg$seed <- 1000 * scenario + s
    
    res <- run_one_replicate(cfg = cfg, demo_stats = demo_stats, real_process = scenario)
    
    meta   <- res$meta
    rep_id <- paste0("S", scenario, "_R", s)
    
    # Truth domains by location
    truth_domains_by_loc <- make_truth_domains_by_loc(
      num_locs = cfg$num_locs,
      t_length = cfg$t_length,
      chosen_locs = meta$chosen_locs,
      real_process = meta$real_process,
      true_idx = 61:80
    )
    
    # --- YOUR METHOD (FocusGP domains) ---
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      method = "DomainSelection_FocusGP_RAW",
      selected_domains_by_loc = res$domains_by_loc_raw,
      truth_domains_by_loc = truth_domains_by_loc,
      meta = meta,
      rep_id = rep_id
    )
    k_domain <- k_domain + 1
    
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      method = "DomainSelection_FocusGP_ISGM_TRUNC",
      selected_domains_by_loc = res$domains_by_loc_isgm_trunc,
      truth_domains_by_loc = truth_domains_by_loc,
      meta = meta,
      rep_id = rep_id
    )
    k_domain <- k_domain + 1
    
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      method = "DomainSelection_FocusGP_ISGM_FULL",
      selected_domains_by_loc = res$domains_by_loc_isgm_full,
      truth_domains_by_loc = truth_domains_by_loc,
      meta = meta,
      rep_id = rep_id
    )
    k_domain <- k_domain + 1
    
    # Baselines (return domain masks too)
    lasso_out <- run_lasso(res$full_dat2, res$sim_scores, cfg$num_locs, cfg$t_length)
    fl_out    <- run_fused_lasso(res$full_dat2, res$sim_scores, cfg$num_locs)
    
    bs13_out  <- run_bspline_lm(res$full_dat2, res$sim_scores, cfg$num_locs, cfg$t_length, bn = 13)
    bs23_out  <- run_bspline_lm(res$full_dat2, res$sim_scores, cfg$num_locs, cfg$t_length, bn = 23)
    selected_locs_FPCA <- run_gp_fpca(
      full_dat2   = res$full_dat2,
      sim_scores  = res$sim_scores,
      num_locs    = cfg$num_locs,
      t_length    = cfg$t_length,
      N           = cfg$nx
    )
    
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      "Lasso", lasso_out$domains_by_loc, truth_domains_by_loc, meta, rep_id
    )
    k_domain <- k_domain + 1
    
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      "FusedLasso", fl_out$domains_by_loc, truth_domains_by_loc, meta, rep_id
    )
    k_domain <- k_domain + 1
    
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      "FLM_BS13", bs13_out$domains_by_loc, truth_domains_by_loc, meta, rep_id
    )
    k_domain <- k_domain + 1
    
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      "FLM_BS23", bs23_out$domains_by_loc, truth_domains_by_loc, meta, rep_id
    )
    k_domain <- k_domain + 1
    
    # FPCA-GP has no domain selection → NA masks per location
    domain_rows[[k_domain]] <- log_domain_method_per_loc(
      "FPCA_GP",
      selected_domains_by_loc = replicate(cfg$num_locs, rep(NA, cfg$t_length), simplify = FALSE),
      truth_domains_by_loc = truth_domains_by_loc,
      meta = meta,
      rep_id = rep_id
    )
    k_domain <- k_domain + 1
    
    # -------------------------
    # Node (location) table
    # -------------------------
    true_locs <- meta$true_locs
    
    node_rows[[k_node]] <- log_method(
      "A1_Full_GP",
      selected_locs = rep(TRUE, cfg$num_locs),
      tau = res$tau_pvals$full,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "A2_DomainOnly",
      selected_locs = rep(TRUE, cfg$num_locs),
      tau = res$tau_pvals$trunc,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "A3_Domain_ISGM",
      selected_locs = res$xi_hat_trunc,
      tau = res$tau_pvals$isgm_trunc,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "A4_ISGM_Full",
      selected_locs = res$xi_hat_full,
      tau = res$tau_pvals$isgm_full,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    # classical node outputs
    node_rows[[k_node]] <- log_method(
      "Lasso",
      selected_locs = lasso_out$selected_locs,
      tau = NA_real_,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "FusedLasso",
      selected_locs = fl_out$selected_locs,
      tau = NA_real_,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "BSpline_13",
      selected_locs = bs13_out$selected_locs,
      tau = NA_real_,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "BSpline_23",
      selected_locs = bs23_out$selected_locs,
      tau = NA_real_,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
    
    node_rows[[k_node]] <- log_method(
      "GP_FPCA",
      selected_locs = selected_locs_FPCA,
      tau = NA_real_,
      truth = true_locs,
      meta = meta,
      rep_id = rep_id
    )
    k_node <- k_node + 1
  }
}
  
domain_results_df <- do.call(rbind, domain_rows)
node_results_df   <- do.call(rbind, node_rows)

print(domain_results_df)
print(node_results_df)

write.csv(domain_results_df, "C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_domain_results_10_LN.csv", row.names = FALSE)
write.csv(node_results_df,   "C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_node_results_10_LN.csv",   row.names = FALSE)


#############################################
#############################################

# =========================================================
# Sensitivity analysis for GP lengthscale (rho)
# =========================================================

rho_grid <- c(1e-2, 0.1, 0.5, 1, 5, 10, 50, 100)
n_reps_per_rho <- 5

sens_rows <- list()
k <- 1

for (rho in rho_grid) {
  
  cfg$lengthscale <- rho   # ← ONLY CHANGE
  
  for (r in seq_len(n_reps_per_rho)) {
    
    cfg$seed <- 10000 + r
    cfg$lengthscale <- rho   # ✅ THIS IS THE ONLY CHANGE YOU NEED
    
    # ---------- NULL ----------
    res_null <- run_one_replicate(
      cfg = cfg,
      demo_stats = demo_stats,
      real_process = 0
    )
    
    sens_rows[[k]] <- data.frame(
      rho = rho,
      replicate = r,
      scenario = "Null",
      p_value = res_null$tau_pvals$isgm_trunc
    )
    k <- k + 1
    
    # ---------- ALTERNATIVE ----------
    res_alt <- run_one_replicate(
      cfg = cfg,
      demo_stats = demo_stats,
      real_process = 1
    )
    
    sens_rows[[k]] <- data.frame(
      rho = rho,
      replicate = r,
      scenario = "Alternative",
      p_value = res_alt$tau_pvals$isgm_trunc
    )
    k <- k + 1
  }
}

# reset for safety
cfg$lengthscale <- NULL

sensitivity_df <- do.call(rbind, sens_rows)
write.csv(sensitivity_df, "C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_rho_sensitivity_isgmtrunc_10_LN.csv", row.names = FALSE)

library(ggplot2)

ggplot(sensitivity_df, aes(x = factor(rho), y = p_value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ scenario, ncol = 1) +
  labs(
    x = "Lengthscale",
    y = "p-value",
    title = "Sensitivity to Lengthscale Choice"
  ) +
  theme_bw()


