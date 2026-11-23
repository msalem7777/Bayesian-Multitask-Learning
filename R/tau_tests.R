# =========================================================
# tau_tests.R
# Tau score test utilities for functional predictors
# =========================================================
# Wraps rho and sigma estimation + tau.score.test
# Requires: hetGP, tau.score.test (from your testing code)

# -------------------------------
# Tau test wrapper
# -------------------------------
#' @export
estimate_sigma2 <- function(GP_df_y) {
  GP_df_y <- as.numeric(GP_df_y)
  N <- length(GP_df_y)
  as.numeric((1 / N) * t(GP_df_y) %*% GP_df_y)
}

#' @export
estimate_rho <- function(GP_df_X) {
  GP_df_X <- as.matrix(GP_df_X)
  rho_vec <- numeric(ncol(GP_df_X))
  
  for (u in seq_len(ncol(GP_df_X))) {
    rho_ph <- darg(NULL, as.matrix(GP_df_X[, u]))
    rho_vec[u] <- (rho_ph$ab[1] - 1) / rho_ph$ab[2]
  }
  
  min(rho_vec, na.rm = TRUE)
}

#' @export
estimate_rho_vec <- function(GP_df_X) {
  GP_df_X <- as.matrix(GP_df_X)
  rho_vec <- numeric(ncol(GP_df_X))
  for (u in seq_len(ncol(GP_df_X))) {
    rho_ph <- darg(NULL, as.matrix(GP_df_X[, u]))
    rho_vec[u] <- (rho_ph$ab[1] - 1) / rho_ph$ab[2]
  }
  rho_vec  # return full vector, not min
}

#' Tau score test with supplied rho vector
#'
#' @return list with statistic and p-value
#' @export
run_tau_test_fixed <- function(Z, y, rho, X = NULL, sigma_vec = NULL){
  
  if(is.null(sigma_vec)){
    sigma_vec <- rep(1, length(y))
  }
  
  pval <- tau.score.test(Z, y, sigma_vec, rho = rho, X = X)
  
  return(pval)
}

#' Tau score test scanning over rho grid
#'
#' @return list with statistic and p-value
#' @export
run_tau_test_scan <- function(Z, y, rho = NULL, X = NULL, sigma_vec = NULL){
  
  if(is.null(sigma_vec)){
    sigma_vec <- rep(1, length(y))
  }
  
  pvals <- tau.score.test2(Z, y, sigma_vec, rho = rho, X = X)
  
  return(pvals)
}

# -------------------------------
# Convenience: build GP matrices
# -------------------------------

#' Stack list of data.frames into matrix
#' @export
stack_locations <- function(dat_list){
  stopifnot(length(dat_list) >= 1)
  
  X <- as.matrix(dat_list[[1]])
  
  if(length(dat_list) > 1){
    for(l in 2:length(dat_list)){
      X <- cbind(X, as.matrix(dat_list[[l]]))
    }
  }
  
  X
}


###############################################################################################
# Score test
#' @export
tau.score.test2 = function(Z, y, sigma_vec, rho = NULL, X=NULL){
  require(hetGP)
  
  if(is.null(rho)){
    rho_vec = c(0.1,0.5,1,1.5,2,3,10,30,50,200)
  } else{
    rho_vec = seq((0.1*min(rho)),(10*max(rho)), length.out = 30)  
  }
  
  Q_vec = P_vec = c()
  
  for(its in 1:length(rho_vec)){
    rho = rho_vec[its]
    S = diag(sigma_vec, nrow = length(y))
    Si = solve(S)
    I = diag(nrow = length(y))
    K = cov_gen(as.matrix(Z), theta=rep(rho, ncol(Z)), type = "Gaussian")
    if((is.null(X))){
      P_0 = Si
      tr.comp = 0.5*sum(diag(P_0%*%K))
      Q.comp = 0.5*t(y)%*%Si%*%K%*%Si%*%(y)
    } else {
      P_0 = Si%*%(I-X%*%solve(t(X)%*%Si%*%X)%*%t(X)%*%Si)
      tr.comp = 0.5*sum(diag(P_0%*%K))
      Q.comp = 0.5*t(y)%*%Si%*%K%*%Si%*%(y)
    }
    
    ITT = sum(diag((P_0%*%K)^2))/2
    ITS = sum(diag(P_0%*%K%*%P_0))/2
    ISS = sum(diag(P_0^2))/2
    ITT_tilde = ITT - ITS%*%solve(ISS)%*%t(ITS)
    
    e_tilde = sum(diag(P_0%*%K))/2
    kappa = ITT/(2*e_tilde)
    df_chisq = (2*e_tilde^2)/ITT_tilde
    
    if((kappa < 0)|(df_chisq<0)){
      kappa = ITT/(2*tr.comp)
      df_chisq = (2*tr.comp^2)/ITT  
    }
    
    score.statistic = Q.comp/kappa
    pval.score = pchisq(score.statistic, df = df_chisq, lower.tail = FALSE)
    
    Q_vec = c(Q_vec, score.statistic)
    P_vec = c(P_vec, pval.score)
  }
  
  # mod_test = summary(lm(P_vec~rho_vec))
  # mod_test = mod_test$fstatistic
  # test.statistic = pf(mod_test[1], mod_test[2], mod_test[3] , lower.tail = FALSE)
  # test.statistic = replace(test.statistic, is.na(test.statistic),1)
  
  # return(1-test.statistic) 
  
  # return(mean(P_vec, na.rm=TRUE))
  
  return(P_vec)
}

#' @export
tau.score.test = function(Z, y, sigma_vec, rho, X=NULL){
  require(hetGP)
  
  S = diag(sigma_vec, nrow = length(y))
  Si = solve(S)
  I = diag(nrow = length(y))
  if(sum((rho>0)*1)==0){
    K = I
  } else {
    K = cov_gen(as.matrix(Z[,as.logical(rho>0)]), theta=rho[as.logical(rho>0)], type = "Gaussian")  
  }
  if((is.null(X))){
    P_0 = Si
    tr.comp = 0.5*sum(diag(P_0%*%K))
    Q.comp = 0.5*t(y)%*%Si%*%K%*%Si%*%(y)
  } else {
    P_0 = Si%*%(I-X%*%solve(t(X)%*%Si%*%X)%*%t(X)%*%Si)
    tr.comp = 0.5*sum(diag(P_0%*%K))
    Q.comp = 0.5*t(y)%*%Si%*%K%*%Si%*%(y)
  }
  
  ITT = sum(diag((P_0%*%K)^2))/2
  ITS = sum(diag(P_0%*%K%*%P_0))/2
  ISS = sum(diag(P_0^2))/2
  ITT_tilde = ITT - ITS%*%solve(ISS)%*%t(ITS)
  
  e_tilde = sum(diag(P_0%*%K))/2
  
  # kappa = as.numeric(ITT_tilde/(2*e_tilde))
  # df_chisq = as.numeric((2*e_tilde^2)/ITT_tilde)
  
  # if((kappa < 0)|(df_chisq<0)){
  #   kappa = ITT/(2*tr.comp)
  #   df_chisq = (2*tr.comp^2)/ITT  
  # }
  
  kappa = ITT/(2*tr.comp)
  df_chisq = (2*tr.comp^2)/ITT  
  
  score.statistic = Q.comp/kappa
  pval.score = pchisq(score.statistic, df = df_chisq, lower.tail = FALSE)
  # pg = seq(1,1000)
  # plot(pg, dchisq(pg, df = df_chisq), type = "l")
  PK <- P_0 %*% K / 2
  PK <- (PK + t(PK)) / 2          # symmetrise to remove floating-point imaginary parts
  eg <- eigen(PK, symmetric = TRUE)
  lb <- Re(eg$values)
  lb <- lb[lb > .Machine$double.eps]   # keep only positive real eigenvalues
  pval.score <- psum.chisq(Q.comp,lb,df=rep(1,length(lb)),nc=rep(0,length(lb)),sigz=0,lower.tail=FALSE,tol=2e-5,nlim=100000,trace=FALSE)
  pval.score <- Re(pval.score)
  
  # If psum.chisq underflows to 0, fall back to the scaled chi-squared approximation
  if (is.na(pval.score) || pval.score <= 0) {
    pval.score <- pchisq(score.statistic, df = df_chisq, lower.tail = FALSE)
  }
  
  return(pval.score)
}