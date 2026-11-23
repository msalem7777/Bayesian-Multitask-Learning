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



