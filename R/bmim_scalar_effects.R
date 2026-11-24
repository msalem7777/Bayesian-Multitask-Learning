# Conjugate scalar effects for the joint BMIM sampler (base R only).
# Under either model beta | noise variance v ~ N(mu, v * B).
# Alternative: v = tau2*g. Null: v = sigma0^2.
# A single shared beta is used for explicit sampling, including b switches.

bmim_scalar_setup <- function(Z, n, beta_prior_mean = NULL,
                              beta_prior_cov = NULL) {
  if (is.null(Z)) {
    if (!is.null(beta_prior_mean) || !is.null(beta_prior_cov))
      stop("Beta prior requires a scalar design Z.")
    return(NULL)
  }
  Z <- as.matrix(Z)
  if (!is.numeric(Z) || nrow(Z) != n || ncol(Z) < 1L ||
      any(!is.finite(Z))) stop("Z must be a finite numeric n-by-S matrix.")
  s <- ncol(Z)
  if (is.null(beta_prior_mean)) beta_prior_mean <- rep(0, s)
  if (is.null(beta_prior_cov)) beta_prior_cov <- diag(100, s)
  mu <- beta_prior_mean
  B <- as.matrix(beta_prior_cov)
  if (!is.numeric(mu) || length(mu) != s || any(!is.finite(mu)) ||
      !is.numeric(B) || !all(dim(B) == c(s,s)) || any(!is.finite(B)) ||
      !isSymmetric(unname(B))) stop("Invalid beta prior mean/covariance.")
  R <- tryCatch(chol(B), error = function(e) NULL)
  if (is.null(R)) stop("Beta prior covariance must be positive definite.")
  nm <- colnames(Z)
  if (is.null(nm)) nm <- paste0("beta", seq_len(s))
  list(Z = Z, s = s, mu = as.numeric(mu), B = B,
       B_inv = chol2inv(R), B_chol = R,
       logdet_B = 2*sum(log(diag(R))),
       B_Z = tcrossprod(Z %*% t(R)), Z_mu = as.numeric(Z %*% mu),
       names = nm)
}

# Return sufficient statistics for integration over the active variance.
# In integrated mode beta has already been averaged out analytically.
# In sampling mode the JOINT density of y and beta must be used, including
# the prior determinant and quadratic form; these depend on g under b=1.
# In both cases the variance is subsequently integrated using the IG prior.
bmim_scalar_terms <- function(y, C, noise_multiplier, scalar, beta = NULL) {
  n <- length(y)
  if (is.null(scalar)) {
    out <- bmim_gp_cholesky_terms(C, y)
    if (!is.null(out)) out$dimension <- n
    return(out)
  }
  if (!is.finite(noise_multiplier) || noise_multiplier <= 0) return(NULL)
  if (is.null(beta)) {
    out <- bmim_gp_cholesky_terms(C + noise_multiplier*scalar$B_Z,
                                  y - scalar$Z_mu)
    if (!is.null(out)) out$dimension <- n
  } else {
    out <- bmim_gp_cholesky_terms(C, y - as.numeric(scalar$Z %*% beta))
    if (!is.null(out)) {
      z <- forwardsolve(t(scalar$B_chol), beta - scalar$mu)
      out$quad <- out$quad + sum(z*z)/noise_multiplier
      out$logdet <- out$logdet + scalar$logdet_B +
        scalar$s*log(noise_multiplier)
      out$dimension <- n + scalar$s
    }
  }
  if (is.null(out) || !is.finite(out$quad) || !is.finite(out$logdet)) return(NULL)
  out
}

# Full normalizing constants are essential for comparing the GP and null.
bmim_scalar_log_evidence <- function(terms, shape, rate) {
  if (is.null(terms)) return(-Inf)
  if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0)
    stop("Variance prior shape and rate must be positive.")
  k <- terms$dimension/2
  lgamma(shape+k)-lgamma(shape)+shape*log(rate)-k*log(2*pi)-
    terms$logdet/2-(shape+k)*log(rate+terms$quad/2)
}

bmim_scalar_alt_terms <- function(y, D, d, g, scalar, beta = NULL) {
  bmim_scalar_terms(y, bmim_gp_covariance(D,d,g), g, scalar, beta)
}

# Under the null C=I, so no n-by-n factorization is needed in sampling
# mode. Integrated mode uses a precomputed result in the MCMC driver.
bmim_scalar_null_terms <- function(y, scalar, beta = NULL) {
  if (is.null(scalar)) return(list(quad=sum(y*y), logdet=0, dimension=length(y)))
  if (is.null(beta)) return(bmim_scalar_terms(y, diag(length(y)), 1, scalar))
  r <- y - as.numeric(scalar$Z %*% beta)
  z <- forwardsolve(t(scalar$B_chol), beta - scalar$mu)
  list(quad=sum(r*r)+sum(z*z), logdet=scalar$logdet_B,
       dimension=length(y)+scalar$s)
}

# Exact Gaussian conditional, using triangular solves rather than C inverse.
# scale is tau2 for b=1 and sigma0^2 for b=0.
bmim_scalar_beta_conditional <- function(y, C, noise_multiplier, scale, scalar) {
  if (is.null(scalar) || !is.finite(scale) || scale <= 0 ||
      !is.finite(noise_multiplier) || noise_multiplier <= 0)
    stop("Invalid scalar conditional inputs.")
  if (is.null(C)) { # Null covariance is identity.
    Xw <- scalar$Z
    yw <- y
  } else {
    R <- chol(C)
    Xw <- forwardsolve(t(R), scalar$Z)
    yw <- forwardsolve(t(R), y)
  }
  precision <- crossprod(Xw) + scalar$B_inv/noise_multiplier
  Rb <- chol(precision)
  rhs <- crossprod(Xw, yw) + scalar$B_inv %*% scalar$mu/noise_multiplier
  mean <- as.numeric(backsolve(Rb, forwardsolve(t(Rb), rhs)))
  list(mean=mean, precision_chol=Rb, scale=scale)
}

bmim_scalar_draw_beta <- function(y, b, D, gp, scales, scalar) {
  post <- if (b == 1L) {
    bmim_scalar_beta_conditional(y, bmim_gp_covariance(D,gp$d,gp$g),
                                 gp$g, scales$tau2, scalar)
  } else {
    bmim_scalar_beta_conditional(y, NULL, 1, scales$sigma2, scalar)
  }
  post$mean + sqrt(post$scale)*as.numeric(backsolve(post$precision_chol,
                                                    rnorm(scalar$s)))
}
