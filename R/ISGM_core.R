# =========================================================
# ISGM_core.R
# Ising Selection via Gibbs / Variational-style Updates
# =========================================================
# This is a self-contained implementation consistent with your usage:
#   ISGM(path_features, W_mat, y, curr_xi, grps_vec, connec, N, VBprop)
#
# Model (conceptual):
#   y ~ f(X * xi) + noise
#   xi ~ Ising prior with adjacency W_mat
#
# The GP parameters are updated by a persistent Gibbs-within-MH kernel.
# The VB association summary supplies an explicit empirical-Bayes
# inclusion-prior factor; deterministic node flips target the full
# conditional GP + Ising + VB product.
# ---------------------------------------------------------


# Load the standalone GP update module from this file's directory.
# This works even when a caller sources ISGM_core.R via an absolute path.
# Inspect active source() frames from innermost to outermost. This also
# works when ISGM_core.R is sourced inside another sourced driver script.
.bmim_source_files <- Filter(
  Negate(is.null),
  lapply(sys.frames(), function(frame) frame$ofile)
)
if (length(.bmim_source_files) == 0L) {
  stop("Source ISGM_core.R from a real file path.")
}
.bmim_current_file <- tail(.bmim_source_files, 1L)[[1L]]
source(file.path(dirname(normalizePath(.bmim_current_file)),
                 "bmim_gp_updates.R"), local = TRUE)
rm(.bmim_current_file, .bmim_source_files)

# -------------------------------
# Helper: group-wise feature aggregation
# -------------------------------

aggregate_by_groups <- function(X, grps_vec){
  # X: N x P design (stacked over locations)
  # grps_vec: list of lists, one per location, giving column groups
  # returns N x L aggregated design

  N <- nrow(X)
  L <- length(grps_vec)

  Z <- matrix(0, N, L)

  col_start <- 1
  for(l in seq_len(L)){
    groups <- grps_vec[[l]]
    cols_l <- unlist(groups)
    Z[,l] <- rowMeans(X[, col_start + cols_l - 1, drop = FALSE])
    col_start <- col_start + max(cols_l)
  }

  Z
}


# -------------------------------
# Helper: logistic function
# -------------------------------
sigmoid = function(x){
  return(1/(1+exp(-x)))
}

# -------------------------------
# Main ISGM function
# -------------------------------

#' Ising Selection with Grouped Features
#' @export
ISGM <- function(path_features, W_mat, y, curr_xi, grps_vec,
                 Z = NULL, N = 1000, connec = 1.0,
                 VBprop = 0, VBsumm = "sq.corr",
                 return_details = FALSE){
  
  # The GP core is base R; infotheo is needed only for mutual information.
  if (VBprop == 1 && VBsumm == "mi" &&
      !requireNamespace("infotheo", quietly = TRUE)) {
    stop("The infotheo package is required for VBsumm=mi.")
  }
  
  num_locs <- length(curr_xi)
  if (num_locs < 1L || any(!curr_xi %in% c(0, 1)) ||
      sum(curr_xi) == 0 || length(N) != 1L || !is.finite(N) ||
      N < 1 || N != as.integer(N)) {
    stop("ISGM requires nonempty binary initial selection and positive integer N.")
  }
  
  # Build mapping of which columns belong to which location (variable truncation lengths)
  map <- build_loc_slices_from_grps(grps_vec)
  loc_slices <- map$loc_slices
  stopifnot(sum(map$lens) == ncol(path_features))
  
  # col_scale[j] = 1 / (number of timepoints in j's location)
  # This implements the 1/(b_p - a_p) normalisation from the kernel spec.
  col_scale <- rep(NA_real_, ncol(path_features))
  for(l in seq_along(loc_slices)){
    n_p <- length(loc_slices[[l]])
    if(n_p > 0) col_scale[loc_slices[[l]]] <- 1 / n_p
  }
  
  # Scale inputs once (always matrix / numeric)
  X_s <- scale(as.matrix(path_features))
  y_s <- as.numeric(scale(y))
  
  # Scalar covariates handling (Fix #5)
  if(is.null(Z)){
    y_tilde <- y_s
    Z_s <- NULL
  } else {
    Z_s <- scale(as.matrix(Z))
    # Minimal, coherent treatment: regress y on Z once and carry residuals through MCMC
    # (If you want full Bayesian sampling of beta, reinstate your commented block properly.)
    fit <- lm.fit(x = cbind(1, Z_s), y = y_s)
    y_tilde <- as.numeric(fit$residuals)
  }
  
  # VB proposal probability p_prop must be per-location (Fix #2)
  if(VBprop == 1){
    if(VBsumm == "sq.corr"){
      g_col <- as.numeric(cor(scale(y), path_features)^2)   # length = ncol(path_features)
    } else if(VBsumm == "mi"){
      g_col <- rep(NA_real_, ncol(path_features))
      for(j in seq_len(ncol(path_features))){
        g_col[j] <- infotheo::mutinformation(infotheo::discretize(y),
                                              infotheo::discretize(path_features[, j]),
                                              method = "emp")
      }
    } else {
      stop("Unknown VBsumm: ", VBsumm)
    }
    
    # Aggregate column-wise g_col to location-wise g_loc (mean over that location's columns)
    g_loc <- sapply(seq_along(loc_slices), function(l){
      idx <- loc_slices[[l]]
      if(length(idx) == 0) return(0)
      mean(g_col[idx], na.rm = TRUE)
    })
    
    g_measure <- data.frame(g_loc)  # num_locs x 1
    
    # IMPORTANT: VBprop() should return a vector of length num_locs
    p_prop <- VBprop(y = g_measure, X = path_features,
                     lambda = 1, h = 0,
                     num_locs = num_locs,
                     corrmat = W_mat,
                     iters = 15, true_dist = 0)
  } else {
    p_prop <- rep(0.5, num_locs)
  }
  
  # Clamp p_prop away from 0/1 to avoid log(0)
  eps <- 1e-8
  p_prop <- pmin(pmax(as.numeric(p_prop), eps), 1 - eps)
  stopifnot(length(p_prop) == num_locs)
  
  # Ising connectivity scaling
  W_mat <- W_mat * connec
  
  # Storage
  cmat <- matrix(NA_real_, nrow = N + 1, ncol = num_locs)
  cmat[1, ] <- curr_xi
  ll_trck <- numeric(0)
  
  # Precompute the same location-normalized distances for GP and node moves.
  D_by_loc <- bmim_gp_precompute_distances(X_s, loc_slices)

  # Explicit, fixed hierarchical prior choices.
  gp_prior <- bmim_gp_default_prior(D_by_loc)
  gp_state <- bmim_gp_initialize(gp_prior)

  for(i in seq_len(N)){

    # ONE persistent GP conditional sweep. Preserve d,g,tau2 across sweeps.
    gp_state <- bmim_gp_update(
      y = y_tilde,
      xi = curr_xi,
      D_by_loc = D_by_loc,
      state = gp_state,
      prior = gp_prior
    )

    # Inner Metropolis flips over locations
    for(j in seq_len(num_locs)){
      
      # Current and proposed xi (flip j)
      prop_xi <- curr_xi
      prop_xi[j] <- 1 - prop_xi[j]
      
      # Optional: forbid selecting nothing
      if(sum(prop_xi) == 0) next
      
      # Log posterior at both node configurations with the SAME GP state.
      # Includes Ising + empirical-Bayes VB + GP determinant + GP quadratic.
      U_current <- bmim_gp_selection_log_target(
        curr_xi, y_tilde, D_by_loc, gp_state, W_mat, p_prop
      )
      U_proposed <- bmim_gp_selection_log_target(
        prop_xi, y_tilde, D_by_loc, gp_state, W_mat, p_prop
      )

      # Deterministic coordinate flips are symmetric: proposal ratio is 1.
      # VB is an explicit factor in the target, not a proposal correction.
      delta_U <- U_proposed - U_current

      prob_flip <- min(1, exp(delta_U))
      flipper <- rbinom(1, 1, prob_flip)
      
      if(flipper == 1){
        curr_xi <- prop_xi
      }
      
      cat("\r", "McMC percentage: ",
          round(((i - 1) * num_locs + j) * 100 / (N * num_locs)),
          "%   ", sep = "")
      flush.console()
    }
    
    cmat[i + 1, ] <- curr_xi
  }
  
  # Return the N post-transition configurations, never the initial state.
  # Callers retain the original matrix interface unless details are requested.
  every_mat <- cmat[seq.int(2L, N + 1L), , drop = FALSE]
  if (!isTRUE(return_details)) return(every_mat)

  # Detailed output makes conditional Bayesian model comparison possible
  # without recomputing the response-dependent VB calculation.
  list(
    samples = every_mat,
    y_tilde = y_tilde,
    D_by_loc = D_by_loc,
    gp_prior = gp_prior,
    W_mat = W_mat,
    p_prop = p_prop,
    gp_last_state = gp_state,
    selected_domains_conditional = TRUE
  )
}

#####################################
# Location selection helper
#####################################
build_loc_slices_from_grps <- function(grps_vec) {
  
  # number of columns contributed by each location after truncation
  lens <- sapply(grps_vec, function(g) length(unname(unlist(g))))
  
  ends   <- cumsum(lens)
  starts <- ends - lens + 1
  
  loc_slices <- lapply(seq_along(lens), function(l) {
    if (lens[l] == 0) integer(0) else starts[l]:ends[l]
  })
  
  list(loc_slices = loc_slices, lens = lens)
}

xi_to_cols <- function(xi, loc_slices) {
  unlist(loc_slices[which(xi == 1)], use.names = FALSE)
}

# Stable quadratic form y' K^{-1} y via Cholesky. Returns Inf if not PD.
quadform_chol <- function(K, y){
  R <- tryCatch(chol(K), error = function(e) NULL)
  if(is.null(R)) return(Inf)
  z <- backsolve(R, y, transpose = TRUE)
  sum(z^2)
}

###############################################################################################
# Selector Function
selector = function(DF, iter_num){
  if(is.na(sum(DF[iter_num,])) & ((iter_num-skipper)==0)) return(DF[1, ])
  if(is.na(sum(DF[iter_num,]))) return(DF[iter_num-skipper, ])
  if(!is.na(sum(DF[iter_num,]))) return(DF[iter_num, ])
}

###############################################################################################
# 0-1 normalization
normalize = function(x, na.rm=TRUE){
  return((x - min(x))/(max(x)-min(x)))
}
###############################################################################################
# Variational Inference
dcorr = function(c2, n, ncp=0){
  if(ncp==0){
    dens = dbeta(c2, 0.5, (n-2)/2)
  } else {
    dens = df(c2, 0.5, (n-2)/2, ncp = ((n*c2)/(1-c2)))
  }
  
  return(dens)
}


VB = function(y, X, corrmat, lambda = 1.0, h=0, iters=50, GP = 0, num_locs, Sigmat = NULL, thresh = 0.9, plot_on=FALSE, true_dist = 0){
  
  ELBO = numeric(iters)
  Hx_mean = numeric(iters)
  
  if(GP == 0){
    library(ClusterR)
    clustr1 = GMM(as.matrix(y), gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10, em_iter = 50)  
    mus = clustr1$centroids
    sds = clustr1$covariance_matrices
    
    emp_mean_num = max(mus)
    emp_sd_num = sqrt(sds[which.max(mus)])
    emp_sd_denom = sqrt(sds[which.min(mus)])
    
    if(true_dist==0){
      logodds = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = emp_sd_num)/(2*dnorm(as.matrix(y), mean = 0, sd = emp_sd_denom)))  
    } else if(true_dist==1){
      # logodds = log(dcorr(as.matrix(y), n = nrow(X), ncp = emp_mean_num)/(dcorr(as.matrix(y), n = nrow(X), ncp = 0)))
      logodds = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = emp_sd_num)/(dcorr(as.matrix(y), n = nrow(X), ncp = 0)))
    }
    
    p1 = sigmoid(logodds)
    mu = 2*p1-1
    
    a = mu + 0.5 * logodds
    
    qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)
    
    if(true_dist==0){
      logp1 = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = sqrt(sds[which.max(mus)])))
      logm1 = log(2*dnorm(as.matrix(y), mean = 0, sd = sqrt(sds[which.min(mus)])))  
    } else if(true_dist==1){
      logp1 = log(dcorr(as.matrix(y), n = nrow(X), ncp = emp_mean_num))
      logm1 = log(dcorr(as.matrix(y), n = nrow(X), ncp = 0))
    }
    
  } else if(GP == 1){
    
    logodds = fX_sig = c()
    fX = list()
    YGP = as.matrix(as.numeric(scale(y)), ncol=1)
    
    for(l in 1:num_locs){
      if((dim(X)[2]/num_locs)>1){
        mod_val = (dim(X)[2]/num_locs)
        XGP = scale(X[,(mod_val*(l-1)+1):(mod_val*(l-1)+mod_val)])
      } else {
        XGP = scale(X[,l])
      }
      
      outg = optim(c(rep(0.1, ncol(XGP)), 0.1*var(YGP)), nlsep, gradnlsep, method="L-BFGS-B", lower=eps, upper=c(rep(3000, ncol(XGP)), 10*var(YGP)), X=XGP, Y=YGP)
      K = plgp::covar.sep(XGP, d=outg$par[1:ncol(XGP)], g=outg$par[ncol(XGP)+1])
      Ki = solve(K)
      tau2hat = drop(t(YGP) %*% Ki %*% YGP / nrow(X))
      KX = covar.sep(XGP, d=outg$par[1:ncol(XGP)], g=0)
      
      fX = c(fX, list(KX %*% Ki %*% YGP))
      fX_sig = c(fX_sig, tau2hat*outg$par[ncol(XGP)+1])
      logodds = c(logodds, (sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log = TRUE))-sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log = TRUE))))
    }
    
    p1 = sigmoid(logodds)
    mu = 2*p1-1
    
    a = mu + 0.5 * logodds
    
    qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)
    
    logp1 = logm1 = c()
    
    for(l in 1:num_locs){
      logp1 = c(logp1, sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log=TRUE)))
      logm1 = c(logm1, sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log=TRUE)))  
    }
  }
  
  #j=1
  for (i in 1:iters){
    # if(i>15){
    #   r = c(1,2,3,4,5,6,7,8,9,10)
    #   model = lm(ELBO[(i-10):(i-1)]~r)
    #   if(abs(coefficients(model)[2])<0.01){
    #     break
    #   }
    # }
    
    muNew = mu
    
    if(GP == 0){
      rlistr = c(sample(seq(1,(length(y))), replace = FALSE))  
    } else if(GP == 1){
      rlistr = c(sample(seq(1,num_locs), replace = FALSE))  
    }
    
    
    for (j in rlistr){
      
      Sbar = sum((corrmat[j,])*muNew)
      if(h==1){
        hp = -sum(muNew)/2*log(sum(muNew)+1)  
        muNew[j] = (1-lambda)*muNew[j] + lambda*tanh(Sbar + 0.5*logodds[j] + hp)
      } else {
        muNew[j] = (1-lambda)*muNew[j] + lambda*tanh(Sbar + 0.5*logodds[j])
      }
      ELBO[i] = ELBO[i] + 0.5*(Sbar * muNew[j])
    }
    
    mu = muNew
    
    a = mu + 0.5 * logodds
    qxp1 = sigmoid(+2*a) #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a) #q_i(x_i=-1)    
    Hx = -qxm1*log(qxm1+1e-10) - qxp1*log(qxp1+1e-10) #entropy        
    
    ELBO[i] = ELBO[i] + sum(qxp1*logp1 + qxm1*logm1) + sum(Hx)
    Hx_mean[i] = mean(Hx)   
    
  }
  
  if(plot_on == TRUE){
    plt_x = seq(1:(i-1))
    plot(plt_x,ELBO[1:(i-1)])
  }
  
  mu_r = round(mu)
  mu_r = replace(mu_r, mu_r == -1, 0)
  
  return(mu_r)
} 

VBprop = function(y, X, corrmat, lambda = 1.0, h=0, iters=50, GP = 0, num_locs, Sigmat = NULL, thresh = 0.9, plot_on=FALSE, true_dist = 0){
  
  ELBO = numeric(iters)
  Hx_mean = numeric(iters)
  
  if(GP == 0){
    library(ClusterR)
    clustr1 = GMM(as.matrix(y), gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10, em_iter = 50)  
    mus = clustr1$centroids
    sds = clustr1$covariance_matrices
    
    emp_mean_num = max(mus)
    emp_sd_num = sqrt(sds[which.max(mus)])
    emp_sd_denom = sqrt(sds[which.min(mus)])
    
    if(true_dist==0){
      logodds = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = emp_sd_num)/(2*dnorm(as.matrix(y), mean = 0, sd = emp_sd_denom)))  
    } else if(true_dist==1){
      # logodds = log(dcorr(as.matrix(y), n = nrow(X), ncp = emp_mean_num)/(dcorr(as.matrix(y), n = nrow(X), ncp = 0)))
      logodds = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = emp_sd_num)/(dcorr(as.matrix(y), n = nrow(X), ncp = 0)))
    }
    
    p1 = sigmoid(logodds)
    mu = 2*p1-1
    
    a = mu + 0.5 * logodds
    
    qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)
    
    if(true_dist==0){
      logp1 = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = sqrt(sds[which.max(mus)])))
      logm1 = log(2*dnorm(as.matrix(y), mean = 0, sd = sqrt(sds[which.min(mus)])))  
    } else if(true_dist==1){
      logp1 = log(dcorr(as.matrix(y), n = nrow(X), ncp = emp_mean_num))
      logm1 = log(dcorr(as.matrix(y), n = nrow(X), ncp = 0))
    }
    
  } else if(GP == 1){
    
    logodds = fX_sig = c()
    fX = list()
    YGP = as.matrix(as.numeric(scale(y)), ncol=1)
    
    for(l in 1:num_locs){
      if((dim(X)[2]/num_locs)>1){
        mod_val = (dim(X)[2]/num_locs)
        XGP = scale(X[,(mod_val*(l-1)+1):(mod_val*(l-1)+mod_val)])
      } else {
        XGP = scale(X[,l])
      }
      
      outg = optim(c(rep(0.1, ncol(XGP)), 0.1*var(YGP)), nlsep, gradnlsep, method="L-BFGS-B", lower=eps, upper=c(rep(3000, ncol(XGP)), 10*var(YGP)), X=XGP, Y=YGP)
      K = plgp::covar.sep(XGP, d=outg$par[1:ncol(XGP)], g=outg$par[ncol(XGP)+1])
      Ki = solve(K)
      tau2hat = drop(t(YGP) %*% Ki %*% YGP / nrow(X))
      KX = covar.sep(XGP, d=outg$par[1:ncol(XGP)], g=0)
      
      fX = c(fX, list(KX %*% Ki %*% YGP))
      fX_sig = c(fX_sig, tau2hat*outg$par[ncol(XGP)+1])
      logodds = c(logodds, (sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log = TRUE))-sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log = TRUE))))
    }
    
    p1 = sigmoid(logodds)
    mu = 2*p1-1
    
    a = mu + 0.5 * logodds
    
    qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)
    
    logp1 = logm1 = c()
    
    for(l in 1:num_locs){
      logp1 = c(logp1, sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log=TRUE)))
      logm1 = c(logm1, sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log=TRUE)))  
    }
  }
  
  #j=1
  for (i in 1:iters){
    # if(i>15){
    #   r = c(1,2,3,4,5,6,7,8,9,10)
    #   model = lm(ELBO[(i-10):(i-1)]~r)
    #   if(abs(coefficients(model)[2])<0.01){
    #     break
    #   }
    # }
    
    muNew = mu
    
    if(GP == 0){
      rlistr = c(sample(seq(1,(length(y))), replace = FALSE))  
    } else if(GP == 1){
      rlistr = c(sample(seq(1,num_locs), replace = FALSE))  
    }
    
    
    for (j in rlistr){
      
      Sbar = sum((corrmat[j,])*muNew)
      if(h==1){
        hp = -sum(muNew)/2*log(sum(muNew)+1)  
        muNew[j] = (1-lambda)*muNew[j] + lambda*tanh(Sbar + 0.5*logodds[j] + hp)
      } else {
        muNew[j] = (1-lambda)*muNew[j] + lambda*tanh(Sbar + 0.5*logodds[j])
      }
      ELBO[i] = ELBO[i] + 0.5*(Sbar * muNew[j])
    }
    
    mu = muNew
    
    a = mu + 0.5 * logodds
    qxp1 = sigmoid(+2*a) #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a) #q_i(x_i=-1)    
    Hx = -qxm1*log(qxm1+1e-10) - qxp1*log(qxp1+1e-10) #entropy        
    
    ELBO[i] = ELBO[i] + sum(qxp1*logp1 + qxm1*logm1) + sum(Hx)
    Hx_mean[i] = mean(Hx)   
    
  }
  
  if(plot_on == TRUE){
    plt_x = seq(1:(i-1))
    plot(plt_x,ELBO[1:(i-1)])
  }
  
  p_r = (mu+1)/2
  
  return(p_r)
} 

###############################################################################################
# Gaussian Kernel Creation
GaussianKernelize = function(path_features, curr_xi, theta){
  
  require(plgp)
  
  holdr = as.matrix(path_features)%*%diag(curr_xi)
  if (sum(colSums(abs(holdr)), na.rm = TRUE)>0){
    holdr = (holdr[, colSums(abs(holdr)) != 0])
  } else {
    holdr = (as.matrix(path_features)[, sample(1:dim(path_features)[2],2)])
  }
  
  kmat = exp(-plgp::distance(holdr)/theta)
  
  return(kmat)
}
