rm(list = ls())
graphics.off()

source("C:/My Files/repos/Bayesian-Multitask-Learning/R/tau_tests.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/run_simulation.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/data_generation.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/domain_selection.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/domain_selection_core.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/correlation_matrix_constructor.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/location_selection.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/ISGM_core.R")
source("C:/My Files/repos/Bayesian-Multitask-Learning/R/evaluation_metrics.R")

# --------------------------------------------------
# Config
# --------------------------------------------------
cfg <- make_default_cfg(
  seed = 1,
  sim_iters = 1,
  sim_setting = "functional_linear",
  nx = 183,
  num_locs = 10,
  num_T = 20,
  chosen_locs = c(2,5,8),
)

demo_stats <- list(
  sim_mean_pat = 110,
  sim_sd_pat   = 10,
  sim_mean_ctr = 100,
  sim_sd_ctr   = 10
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
run_bspline_lm <- function(
    full_dat2,
    sim_scores,
    num_locs,
    t_length,
    bn
) {
  library(fda)
  
  N <- length(sim_scores)
  
  # --- build dp_mat exactly as original ---
  dp_mat <- matrix(0, nrow = N)
  basisobj <- create.bspline.basis(c(1, t_length), bn)
  
  for (l in seq_len(num_locs)) {
    fd_dat <- t(full_dat2[[l]])
    colnames(fd_dat) <- seq_len(N)
    rownames(fd_dat) <- seq_len(t_length)
    
    b_fd <- smooth.basisPar(
      seq_len(t_length),
      fd_dat,
      basisobj
    )$fd
    
    dp_mat <- cbind(dp_mat, t(as.matrix(b_fd$coefs)))
  }
  
  dp_mat <- dp_mat[, 2:ncol(dp_mat)]
  colnames(dp_mat) <- paste0("X", seq_len(ncol(dp_mat)))
  
  # --- linear model ---
  fit <- lm(scale(sim_scores) ~ dp_mat)
  fit_summ <- summary(fit)
  
  pvals <- fit_summ$coefficients[
    2:nrow(fit_summ$coefficients),
    4
  ]
  
  # --- NA handling EXACTLY like original ---
  pvals_nan <- seq_along(pvals) * is.na(pvals)
  pvals[pvals_nan] <- 1
  
  if (length(pvals) != ncol(dp_mat)) {
    pvals <- c(pvals, rep(1, ncol(dp_mat) - length(pvals)))
  }
  
  # --- location-level selection (bn blocks) ---
  selected_locs <- logical(num_locs)
  
  for (lnum in seq_len(num_locs)) {
    idx <- ((lnum - 1) * bn + 1):(lnum * bn)
    selected_locs[lnum] <- min(pvals[idx]) < 0.05
  }
  
  selected_locs
}


# --------------------------------------------------
# Fused Lasso
# --------------------------------------------------
run_fused_lasso <- function(full_dat2, sim_scores, num_locs,
                            lambda1 = 0.1, lambda2 = 1, thresh = 0.005) {
  
  library(extlasso)
  
  selected_locs <- logical(num_locs)
  
  for (l in seq_len(num_locs)) {
    fit <- extlasso::fusedlasso(as.matrix(full_dat2[[l]]),
                      scale(sim_scores),
                      lambda1 = lambda1,
                      lambda2 = lambda2,
                      maxiter = 2e5)
    
    selected_locs[l] <- any(abs(fit$coef) > thresh)
  }
  
  selected_locs
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
  
  selected_locs <- logical(num_locs)
  for (l in seq_len(num_locs)) {
    idx <- ((l - 1) * t_length + 1):(l * t_length)
    selected_locs[l] <- any(beta[idx] != 0)
  }
  
  selected_locs
}

# --------------------------------------------------
# Storage
# --------------------------------------------------
rows <- list()

# --------------------------------------------------
# Loop over replicates
# --------------------------------------------------
for (s in seq_len(cfg$sim_iters)) {
  
  cfg$seed <- s
  
  res <- run_one_replicate(
    cfg = cfg,
    demo_stats = demo_stats,
    real_process = 1
  )
  
  true_locs <- res$meta$true_locs
  meta      <- res$meta
  rep_id    <- s
  
  selected_locs_FPCA <- run_gp_fpca(
    full_dat2 = res$full_dat2,
    sim_scores = res$sim_scores,
    num_locs = cfg$num_locs,
    t_length = t_length,
    N = cfg$nx
  )
  
  selected_locs_B20 <- run_bspline_lm(
    full_dat2   = full_dat2,
    sim_scores  = sim_scores,
    num_locs    = cfg$num_locs,
    t_length    = 100,
    bn          = 23
  )
  
  selected_locs_B13 <- run_bspline_lm(
    full_dat2   = full_dat2,
    sim_scores  = sim_scores,
    num_locs    = cfg$num_locs,
    t_length    = 100,
    bn          = 13
  )
  
  selected_locs_FL <- run_fused_lasso(
    full_dat2, sim_scores, cfg$num_locs
  )
  
  selected_locs_CL <- run_lasso(
    full_dat2, sim_scores, cfg$num_locs, t_length
  )
  
  # -------------------------------
  # Method A: GP-based comparisons
  # -------------------------------
  
  rows[[length(rows)+1]] <- log_method(
    method = "A1_Full_GP",
    selected_locs = rep(TRUE, cfg$num_locs),
    tau = res$tau_pvals$full,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    method = "A2_DomainOnly",
    selected_locs = rep(TRUE, cfg$num_locs),
    tau = res$tau_pvals$trunc,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    method = "A3_Domain_ISGM",
    selected_locs = res$xi_hat_trunc,
    tau = res$tau_pvals$isgm_trunc,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    method = "A4_ISGM_Full",
    selected_locs = res$xi_hat_full,
    tau = res$tau_pvals$isgm_full,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  # ------------------------------------------------
  # Classical baselines (placeholders for now)
  # ------------------------------------------------
  # These can be filled once you re-enable the blocks
  
  rows[[length(rows)+1]] <- log_method(
    "GP_FPCA",
    selected_locs = rep(NA, cfg$num_locs),
    tau = NA_real_,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    "BSpline_20",
    selected_locs = rep(NA, cfg$num_locs),
    tau = NA_real_,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    "BSpline_40",
    selected_locs = rep(NA, cfg$num_locs),
    tau = NA_real_,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    "FusedLasso",
    selected_locs = rep(NA, cfg$num_locs),
    tau = NA_real_,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
  
  rows[[length(rows)+1]] <- log_method(
    "Lasso",
    selected_locs = rep(NA, cfg$num_locs),
    tau = NA_real_,
    truth = true_locs,
    meta = meta,
    rep_id = rep_id
  )
}

# --------------------------------------------------
# Bind + save
# --------------------------------------------------
results_df <- do.call(rbind, rows)

write.csv(
  results_df,
  "simulation_comparison_results.csv",
  row.names = FALSE
)

# quick check
print(results_df)



############################################
cfg <- list(
  seed = 1,
  nx = 183,
  num_locs = 10,
  num_T = 20,
  gt_sigma_pop = 0.25,
  gt_tau = 4,
  gt_sigma_u = 4,
  chosen_locs = c(2,5,8),
  chosen_one = 1,
  unchosen_one = 2,
  sim_setting = "real_or_noise"
)

# IMPORTANT: t_length must match block construction:
# prepend (20) + cycle1 (20) + append1 (20) + cycle2 (20) + append2 (20) = 100
cfg$t_length <- 5 * cfg$num_T

str(cfg)

cfg$sim_setting <- "real_or_noise"

out1 <- simulate_one_dataset(
  cfg,
  sim_mean_pat = 110, sim_sd_pat = 10,
  sim_mean_ctr = 100, sim_sd_ctr = 10,
  real_process = 1
)

str(out1, max.level = 1)

lapply(out1$full_dat2, dim)

#######
loc_sig <- 2
loc_null <- 1

mat_sig <- as.matrix(out1$full_dat2[[loc_sig]])
mat_null <- as.matrix(out1$full_dat2[[loc_null]])

matplot(t(mat_sig[1:30, ]), type = "l", lty = 1,
        main = paste("Signal location", loc_sig),
        xlab = "Time", ylab = "Value")

matplot(t(mat_null[1:30, ]), type = "l", lty = 1,
        main = paste("Null location", loc_null),
        xlab = "Time", ylab = "Value")

######
case_idx <- which(out1$true_status == 1)
ctrl_idx <- which(out1$true_status == 0)

matplot(t(mat_sig[case_idx[1:15], ]), type = "l", col = "red",
        main = "Signal location – CASES", xlab = "Time", ylab = "Value")

matplot(t(mat_sig[ctrl_idx[1:15], ]), type = "l", col = "blue",
        main = "Signal location – CONTROLS", xlab = "Time", ylab = "Value")

#########
cfg$sim_setting <- "functional_linear"

out2 <- simulate_one_dataset(
  cfg,
  sim_mean_pat = 110, sim_sd_pat = 10,
  sim_mean_ctr = 100, sim_sd_ctr = 10,
  real_process = 1
)

summary(out2$sim_scores)

##########
out2 <- simulate_one_dataset(
  cfg,
  sim_mean_pat = 110, sim_sd_pat = 10,
  sim_mean_ctr = 100, sim_sd_ctr = 10,
  real_process = 1
)

summary(out2$sim_scores)

############
signal_mat <- Reduce("+", lapply(out2$full_dat2[out2$chosen_locs], as.matrix))

signal_window <- rowMeans(signal_mat[, 61:80])

cor(signal_window, out2$sim_scores)

#############
plot(signal_window, out2$sim_scores,
     xlab = "Mean functional signal (true window)",
     ylab = "Outcome",
     main = "Functional-linear coupling check")
abline(lm(out2$sim_scores ~ signal_window), col = "red")

#########
X_mat <- stack_locations(out2$full_dat2)
dim(X_mat)

###########
pvals_scan <- run_tau_test_scan(
  Z = X_mat,
  y = out2$sim_scores
)

str(pvals_scan)
summary(pvals_scan)


plot(pvals_scan, type="b", xlab="rho grid index", ylab="p-value",
     main="Tau score test p-values across rho grid")
abline(h=0.05, lty=2, col="red")


out_null <- simulate_one_dataset(
  cfg,
  sim_mean_pat = 110, sim_sd_pat = 10,
  sim_mean_ctr = 100, sim_sd_ctr = 10,
  real_process = 0
)

X_null <- stack_locations(out_null$full_dat2)

pvals_null <- run_tau_test_scan(
  Z = X_null,
  y = out_null$sim_scores
)

summary(pvals_null)

plot(pvals_null, type="b", xlab="rho grid index", ylab="p-value",
     main="Null tau test p-values across rho grid")
abline(h=0.05, lty=2, col="red")


summary(pvals_scan)
summary(pvals_null)

#######
dom <- run_domain_selection(out2$full_dat2, out2$sim_scores, cfg$t_length, cfg$num_locs)
sapply(dom$final_data, ncol)
lengths(dom$grps_vec)

W_real <- build_W_matrix(
  mode = "data",
  full_dat2 = out2$full_dat2,
  t_length = cfg$t_length
)

image(W_real$W_func, main = "Estimated functional correlation (real-data mode)")
hist(W_real$W_func[upper.tri(W_real$W_func)], breaks = 30)

W_sim <- build_W_matrix(
  mode = "simulation",
  full_dat2 = out2$full_dat2,
  chosen_locs = out2$chosen_locs,
  t_length = cfg$t_length
)

W_sim$W_func
W_sim$W_func_full

image(W_sim$W_func, main = "W_func (boosted chosen_locs)")
image(W_sim$W_func_full, main = "W_func_full (weaker)")


dom <- run_domain_selection(
  full_dat2 = out2$full_dat2,
  sim_scores = out2$sim_scores,
  t_length = cfg$t_length,
  num_locs = cfg$num_locs
)

sapply(dom$final_data, ncol)

# ----------------------------
# 1) Build truncated design EXACTLY like original (but numeric-safe)
# ----------------------------
full_dat2 <- out2$full_dat2  # keep your original naming style

num_locs <- length(dom$grps_vec)

final_data <- vector("list", num_locs)
for (l in seq_len(num_locs)) {
  df_trunc <- as.data.frame(full_dat2[[l]])
  df_trunc <- df_trunc[, unname(unlist(dom$grps_vec[[l]])), drop = FALSE]
  final_data[[l]] <- df_trunc
}

# IMPORTANT: make a numeric, wide design matrix/data.frame
# (as.data.frame(final_data) can create list-columns -> breaks colMeans/scale)
X_trunc <- as.data.frame(do.call(cbind, final_data))

grps_trunc <- dom$grps_vec

# Sanity check: ISGM_core mapping must match columns
map_tr <- build_loc_slices_from_grps(grps_trunc)
cat("sum(lens) =", sum(map_tr$lens), " ncol(X_trunc) =", ncol(X_trunc), "\n")
stopifnot(sum(map_tr$lens) == ncol(X_trunc))


# ----------------------------
# 2) Build NON-truncated design (full) + matching groups
# ----------------------------
final_data_full <- lapply(full_dat2, function(df) {
  as.data.frame(df)[, seq_len(cfg$t_length), drop = FALSE]
})

X_full <- as.data.frame(do.call(cbind, final_data_full))

# Make grps_full in the SAME FORMAT as grps_trunc (list of groups per location)
grps_full <- lapply(final_data_full, function(df) {
  list(seq_len(ncol(df)))
})

map_fl <- build_loc_slices_from_grps(grps_full)
cat("sum(lens) =", sum(map_fl$lens), " ncol(X_full) =", ncol(X_full), "\n")
stopifnot(sum(map_fl$lens) == ncol(X_full))


# ----------------------------
# 3) Run ISGM twice (trunc + full)
# ----------------------------
set.seed(1)
xi0 <- rep(1, cfg$num_locs)

xi_chain_trunc <- ISGM(
  path_features = X_trunc,
  W_mat = W_sim$W_func,
  y = out2$sim_scores,
  curr_xi = xi0,
  grps_vec = grps_trunc,
  connec = 1,
  N = 300,
  VBprop = 1
)

xi_chain_full <- ISGM(
  path_features = X_full,
  W_mat = W_sim$W_func_full,
  y = out2$sim_scores,
  curr_xi = xi0,
  grps_vec = grps_full,
  connec = 1,
  N = 300,
  VBprop = 1
)

dim(xi_chain_trunc)
dim(xi_chain_full)

# quick look at last-iter selection
tail(xi_chain_trunc, 1)
tail(xi_chain_full, 1)


# ============================================================
# TAU TEST RECONSTRUCTION (Full vs Domain vs Domain+ISGM)
# Drop this block AFTER you have:
#   - out2
#   - dom (from run_domain_selection)
#   - xi_chain_trunc (from ISGM on truncated)
#   - X_full (full design, numeric-safe)
#   - normalize(), darg(), tau.score.test() in scope
# ============================================================

# ----------------------------
# 0) Objects + terminology alignment (match original naming)
# ----------------------------
full_dat2 <- out2$full_dat2
N <- length(out2$sim_scores)

# IMPORTANT: original code uses final_int_df as the (domain-truncated) wide design
# In your numeric-safe refactor, that is exactly X_trunc
# If you haven't already set it, do it here:
if (!exists("X_trunc")) {
  stop("X_trunc not found. Build X_trunc with do.call(cbind, final_data) before running tau block.")
}
final_int_df <- X_trunc

# ----------------------------
# 1) Intercept-only residualization (exact original)
# ----------------------------
X_s <- matrix(1, nrow = N, ncol = 1)
modl_s <- lm(out2$sim_scores ~ X_s - 1)

GP_df_y <- as.numeric(out2$sim_scores - modl_s$fitted.values)
sigma_vec <- as.numeric((1 / N) * t(GP_df_y) %*% GP_df_y)

# ----------------------------
# 2) rho estimator (exact original)
# ----------------------------
estimate_rho <- function(GP_df_X) {
  GP_df_X <- as.matrix(GP_df_X)
  rho_vec <- numeric(ncol(GP_df_X))
  
  for (u in seq_len(ncol(GP_df_X))) {
    rho_ph <- darg(NULL, as.matrix(GP_df_X[, u]))
    rho_vec[u] <- (rho_ph$ab[1] - 1) / rho_ph$ab[2]
  }
  
  min(rho_vec, na.rm = TRUE)
}

# ----------------------------
# 3) TAU #1: Full (no selection)
#     Original used: as.data.frame(full_dat2) then normalize
#     Your equivalent is X_full (numeric wide design)
# ----------------------------
if (!exists("X_full")) {
  stop("X_full not found. Build X_full with do.call(cbind, final_data_full) before running tau block.")
}

GP_df_X_full <- apply(as.matrix(X_full), 2, normalize)
rho_full <- estimate_rho(GP_df_X_full)

tau_full <- tau.score.test(
  GP_df_X_full,
  GP_df_y,
  sigma_vec,
  rho = rho_full,
  X = X_s
)

# ----------------------------
# 4) TAU #2: Domain selection only ("Integral")
#     Original used: apply(final_int_df, 2, normalize)
# ----------------------------
GP_df_X_int <- apply(as.matrix(final_int_df), 2, normalize)
rho_int <- estimate_rho(GP_df_X_int)

tau_int <- tau.score.test(
  GP_df_X_int,
  GP_df_y,
  sigma_vec,
  rho = rho_int,
  X = X_s
)

# ----------------------------
# 5) TAU #3: Domain + ISGM (Integral + Ising selection)
#     Must mask by LOCATION blocks, not raw columns.
#     Use last sample of xi_chain_trunc like original "selection = curr_xi_dfs[4,]"
# ----------------------------
if (!exists("xi_chain_trunc")) {
  stop("xi_chain_trunc not found. Run ISGM on truncated data before running tau block.")
}

selection_trunc <- as.logical(tail(xi_chain_trunc, 1))  # length = num_locs
num_locs <- length(selection_trunc)

# Build location-block map for final_int_df:
# Each location contributes: length(unlist(dom$grps_vec[[l]])) columns
block_sizes <- sapply(dom$grps_vec, function(g) length(unlist(g)))
stopifnot(length(block_sizes) == num_locs)

block_idx <- rep(seq_len(num_locs), times = block_sizes)
stopifnot(length(block_idx) == ncol(final_int_df))

if (sum(selection_trunc) > 0) {
  GP_df_X_ising <- final_int_df[, block_idx %in% which(selection_trunc), drop = FALSE]
  GP_df_X_ising <- apply(as.matrix(GP_df_X_ising), 2, normalize)
  rho_ising <- estimate_rho(GP_df_X_ising)
} else {
  # match original fallback behavior
  GP_df_X_ising <- matrix(0, nrow = N, ncol = 1)
  rho_ising <- 0.1
}

tau_ising <- tau.score.test(
  GP_df_X_ising,
  GP_df_y,
  sigma_vec,
  rho = rho_ising,
  X = X_s
)

# ----------------------------
# 6) Report
# ----------------------------
print(list(
  tau_full  = tau_full,
  tau_int   = tau_int,
  tau_ising = tau_ising,
  rho_full  = rho_full,
  rho_int   = rho_int,
  rho_ising = rho_ising,
  selection_trunc = selection_trunc
))

print(data.frame(
  Method = c("Full (no selection)", "Domain selection only", "Domain + ISGM"),
  Tau = c(tau_full, tau_int, tau_ising),
  Rho = c(rho_full, rho_int, rho_ising),
  Selected_Locs = c(
    paste0("ALL (", num_locs, ")"),
    paste0("ALL (", num_locs, ")"),
    paste(which(selection_trunc), collapse = ",")
  )
))
