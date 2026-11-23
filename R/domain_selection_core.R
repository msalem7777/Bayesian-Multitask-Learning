# =========================================================
# domain_selection_core.R
# Core GP optimization utilities used by FocusGP
# =========================================================

# This file contains:
#   - squared exponential kernel
#   - GP marginal likelihood
#   - GP prediction
#   - Expected Improvement (EI)
#   - k-fold CV objective
#   - FocusGP interval search wrapper
#
# Dependencies: MASS, mvtnorm

# ------------------------------- 
# k-fold cross validation loss 
# -------------------------------

FCV = function(X, Z, a, s, nfold, lambda = 0.1){
  # initialize our covariates
  X = as.data.frame(X)
  
  a = max(c(a,1))
  a = min(c(a,90))
  s = max(c(1,s))
  # s = min(c(s, (100-a)))
  
  # X1 = apply(X[,a:min(c(100,(a+s)))],1, mean, na.rm=TRUE)
  # X2 = apply(X[,a:min(c(100,(a+s)))],1, var, na.rm=TRUE)
  # X_df = data.frame(X1=X1, X2=X2)
  X_df = X[,a:min(c(100,(a+s)))]
  # fit best GP given a and s and compute average validation error
  rmses = c()
  for(k in 1:nfold){
    ind_s = (k-1)*round(nrow(X_df)/nfold)+1
    ind_e = round(nrow(X_df)/nfold*k)
    XT = X_df[ind_s:ind_e,]
    YT = Z[ind_s:ind_e]
    XV = X_df[-c(ind_s:ind_e),]
    YV = Z[-c(ind_s:ind_e)]
    gpiFCV <- newGPsep(XT, YT, d=0.1, g=0.1*var(YT), dK=TRUE)
    mle <- mleGPsep(gpiFCV, param="d", tmin=0.1, tmax=c(3000))
    p <- predGPsep(gpiFCV, XV)
    deleteGPsep(gpiFCV)
    rmses = c(rmses, (sqrt(mean((YV - p$mean)^2)) + lambda*((s-20)^2))) + 99999*(s>(100-a))
  }
  
  ValErr = mean(rmses, na.rm=TRUE)
}

# -------------------------------
# Expected Improvement
# -------------------------------

eps <- sqrt(.Machine$double.eps)

EI = function(gpi, x, fmin, pred=predGPsep){
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  p <- pred(gpi, x, lite=TRUE)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d/sigma
  ei <- d*pnorm(dn) + sigma*dnorm(dn)
  return(ei)
}

obj.EI = function(x, fmin, gpi, pred=predGPsep){
  - EI(gpi, x, fmin, pred)
}

EI.search = function(X, y, gpi, pred=predGPsep, multi.start=1, tol=eps){
  m <- which.min(y)
  fmin <- y[m]
  start <- matrix(X[m,], nrow=1)
  if(multi.start > 1)
    start <- rbind(start, randomLHS(multi.start - 1, ncol(X)))
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  for(i in 1:nrow(start)) {
    if(EI(gpi, start[i,], fmin) <= tol) { out <- list(value=-Inf); next }
    out <- optim(start[i,], obj.EI, method="L-BFGS-B",lower=c(1,5), upper = c(95, 50), gpi=gpi, pred=pred, fmin=fmin)
    xnew[i,] <- c(out$par, -out$value)
  }
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c("s1", "s2", "x1", "x2", "val")
  solns <- solns[solns$val > tol,]
  return(solns)
}


# -------------------------------
# FocusGP wrapper
# -------------------------------

#' Bayesian optimization over time intervals
#' @export
FocusGP = function(X, Z, nfold, lambda, ninit=12, end = 30){
  
  library(lhs)
  library(laGP)
  
  A_init = randomLHS(ninit, 1)
  A_init = round(A_init*(ncol(X)-11))+1
  S_init = round(randomLHS(ninit, 1)*(ncol(X)-(A_init+10))) + 10
  P_init = cbind(A_init, S_init)
  
  # Randomizing data order for validation splits
  rand = sample(nrow(X))
  X = X[rand, ]
  Z = Z[rand]
  
  ValErr = c()
  for(j in 1:nrow(P_init)){
    ValErr = c(ValErr, FCV(X, Z, A_init[j], S_init[j], lambda = 0.1, nfold = nfold))
  }
  
  # Build GP on domain variables
  gpi <- newGPsep(P_init, ValErr, d=0.1, g=1e-6, dK=TRUE)
  da <- darg(list(mle=TRUE, max=0.5), randomLHS(1000, 1))
  da$max = 9999999999999
  mleGPsep(gpi, param="d", tmin=da$min, tmax=(da$max), ab=da$ab)
  # propagate validation error to acquisition function and update a and s
  maxei <- c()
  for(i in (ninit+1):end) {
    solns <- EI.search(P_init, ValErr, gpi)
    m <- which.max(solns$val)
    maxei <- c(maxei, solns$val[m])
    xnew <- as.matrix(solns[m,3:4])
    ynew <- FCV(X, Z, xnew[1], xnew[2], lambda = 0.1, nfold = nfold)
    updateGPsep(gpi, xnew, ynew)
    mleGPsep(gpi, param="d", tmin=da$min, tmax=(da$max), ab=da$ab)
    P_init <- unname(rbind(P_init, xnew))
    ValErr <- c(ValErr, ynew)
  }
  ## clean up and return
  deleteGPsep(gpi)
  
  return(list(X=P_init, y=ValErr, maxei=maxei))
}

nlsep <- function(par, X, Y){
  theta <- par[1:ncol(X)]
  g <- par[ncol(X)+1]
  n <- length(Y)
  K <- covar.sep(X, d=theta, g=g)
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  return(-ll)
}

gradnlsep <- function(par, X, Y){
  theta <- par[1:ncol(X)]
  g <- par[ncol(X)+1]
  n <- length(Y)
  K <- covar.sep(X, d=theta, g=g)
  Ki <- solve(K)
  KiY <- Ki %*% Y
  ## loop over theta components
  dlltheta <- rep(NA, length(theta))
  for(k in 1:length(dlltheta)) {
    dotK <- K * distance(X[,k])/(theta[k]^2)
    dlltheta[k] <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) -
      (1/2)*sum(diag(Ki %*% dotK))
  }
  ## for g
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  return(-c(dlltheta, dllg))
}
