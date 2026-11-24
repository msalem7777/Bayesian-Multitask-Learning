# Scalar-effect tests. Base R only. Run from R/:
# Rscript -e "source('checks/bmim_scalar_effects_smoke.R')"
files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
root <- dirname(dirname(normalizePath(tail(files,1L)[[1L]])))
source(file.path(root,"bmim_gp_updates.R"))
source(file.path(root,"bmim_global_gp_test.R"))
source(file.path(root,"bmim_joint_domains.R"))

# Independent multivariate Student-t log density obtained from a Gaussian
# with covariance v*C and v~IG(a,r). Do not call the new likelihood helpers.
ref_log <- function(y, C, a, r) {
  n <- length(y)
  logdet <- as.numeric(determinant(C, logarithm=TRUE)$modulus)
  q <- as.numeric(crossprod(y, solve(C,y)))
  lgamma(a+n/2)-lgamma(a)+a*log(r)-n/2*log(2*pi)-logdet/2-
    (a+n/2)*log(r+q/2)
}
set.seed(8301)
n <- 7L
z <- seq(-1,1,length.out=n)
X <- scale(cbind(z+0.15*sin(4*z), z^2, sin(3*z), cos(4*z)))
y <- 0.7+0.9*z+0.3*sin(3*z)+rnorm(n,sd=.25)
Z <- cbind(intercept=1,z=z)
mu <- c(.1,-.2)
B <- matrix(c(1.4,.3,.3,.8),2L)
s <- bmim_scalar_setup(Z,n,mu,B)
stopifnot(isTRUE(all.equal(s$B_Z,Z%*%B%*%t(Z),tolerance=1e-12)))
loc <- list(1:2,3:4)
D <- bmim_gp_precompute_distances(X,loc)[[1L]]
d <- .7; g <- .4; a <- 2.3; r <- .8
C <- exp(-D/d)+g*diag(n)
beta <- c(.6,.4)
# Verify full constants for both active models, unequal variance priors,
# nonzero beta means and non-diagonal coefficient covariance.
for (v in list(list(C=C,m=g),list(C=diag(n),m=1))) {
  C0 <- v$C; m <- v$m
  expected <- ref_log(y-Z%*%mu,C0+m*Z%*%B%*%t(Z),a,r)
  got <- bmim_scalar_log_evidence(bmim_scalar_terms(y,C0,m,s),a,r)
  stopifnot(abs(expected-got)<1e-10)
  # Joint Gaussian for (y,beta) before integrating its common scale.
  joint <- rbind(cbind(C0+m*Z%*%B%*%t(Z),m*Z%*%B),
                 cbind(m*B%*%t(Z),m*B))
  expected_joint <- ref_log(c(y-Z%*%mu,beta-mu),joint,a,r)
  got_joint <- bmim_scalar_log_evidence(bmim_scalar_terms(y,C0,m,s,beta),a,r)
  stopifnot(abs(expected_joint-got_joint)<1e-10)
  if (m==1) {
    stopifnot(abs(bmim_scalar_log_evidence(bmim_scalar_null_terms(y,s),a,r)-got)<1e-10,
              abs(bmim_scalar_log_evidence(bmim_scalar_null_terms(y,s,beta),a,r)-got_joint)<1e-10)
  }
  # Independent conditional Gaussian from standard block conditioning.
  post <- bmim_scalar_beta_conditional(y,C0,m,.9,s)
  cross <- m*B%*%t(Z)
  marg <- C0+m*Z%*%B%*%t(Z)
  ref_mean <- as.numeric(mu+cross%*%solve(marg,y-Z%*%mu))
  ref_cov <- .9*(m*B-cross%*%solve(marg,t(cross)))
  stopifnot(max(abs(post$mean-ref_mean))<1e-10,
            max(abs(.9*chol2inv(post$precision_chol)-ref_cov))<1e-10)
}
cat("SCALAR MARGINAL / JOINT LIKELIHOOD AND GAUSSIAN CONDITIONAL: PASS\n")

# Exact posterior for fixed d,g, all three nonempty node configurations and
# all nine interval pairs (length=1 or 2), using an independently built
# normalized node prior and interval prior. Null and alternative IG differ.
W <- matrix(c(0,.15,.15,0),2L)
p_prop <- c(.58,.44); pi_b <- .45
prior <- bmim_gp_default_prior(bmim_gp_precompute_distances(X,loc),
  alpha_tau=4.6,q_tau=1.6,null_sigma_shape=1.8,null_sigma_rate=1.2)
xis <- rbind(c(1,0),c(0,1),c(1,1))
w <- apply(xis,1,function(x) exp(.5*sum((2*x-1)*(W%*%(2*x-1))))*
  prod(p_prop^x*(1-p_prop)^(1-x)))
w <- w/sum(w)
intervals <- list(1L,2L,1:2); mass <- c(.25,.25,.5)
log0 <- ref_log(y-Z%*%mu,diag(n)+Z%*%B%*%t(Z),
                prior$null_sigma_shape,prior$null_sigma_rate)
rows <- list()
for (j in 1:3) for (k in 1:3) for (l in 1:3) {
  idx <- list(intervals[[k]],intervals[[l]])
  D_ref <- matrix(0,n,n)
  for (p in 1:2) if (xis[j,p]==1) {
    for (t in loc[[p]][idx[[p]]])
      D_ref <- D_ref+outer(X[,t],X[,t],"-")^2/length(idx[[p]])
  }
  log1 <- ref_log(y-Z%*%mu,exp(-D_ref/d)+g*(diag(n)+Z%*%B%*%t(Z)),
                  prior$tau_shape,prior$tau_rate)
  rows[[length(rows)+1L]] <- c(xis[j,],
    (1-pi_b)*w[j]*mass[k]*mass[l]*exp(log0),
    pi_b*w[j]*mass[k]*mass[l]*exp(log1))
}
ref <- do.call(rbind,rows)
ref_gp <- sum(ref[,4])/sum(ref[,3:4])
ref_nodes <- colSums(ref[,1:2]*ref[,4])/sum(ref[,4])
args <- list(X_s=X,y_tilde=y,loc_slices=loc,W_mat=W,p_prop=p_prop,
  max_domain_length=2L,n_iter=6000L,burnin=1000L,n_chains=4L,
  pi_b=pi_b,gp_prior=prior,fixed_gp=list(d=d,g=g),Z=Z,
  beta_prior_mean=mu,beta_prior_cov=B,seed=8302L)
for (mode in c("integrated","sampled")) {
  fit <- do.call(bmim_joint_domain_mcmc,c(args,list(beta_mode=mode,save_beta=TRUE)))
  stopifnot(fit$beta_mode==mode,fit$beta_draws_saved,
            abs(fit$posterior_gp_included-ref_gp)<.055,
            max(abs(fit$posterior_node_included_given_gp-ref_nodes))<.07,
            all(vapply(fit$chains,function(ch) all(is.finite(ch$beta)),logical(1))))
  cat(sprintf("SCALAR %s FIXED-GP ENUMERATION: PASS (GP %.4f vs %.4f)\n",
              mode,fit$posterior_gp_included,ref_gp))
}
# Exercise both modes with real d/g updates and the feature wrapper,
# with and without global b. Slice and random-walk g are both covered.
for (mode in c("integrated","sampled")) for (global in c(FALSE,TRUE)) {
  pr <- prior; pr$g_update_method <- if (global) "slice" else "rw"
  fit <- bmim_joint_domains_from_features(X,W,y,loc,p_prop,2L,
    Z=matrix(z,ncol=1),beta_mode=mode,save_beta=TRUE,
    n_iter=60L,burnin=10L,n_chains=2L,pi_b=if (global) pi_b else NULL,
    gp_prior=pr,seed=8303L)
  stopifnot(fit$beta_mode==mode,fit$response_preprocessing$scale==1,
            identical(fit$beta_design_names,c("(Intercept)","Z1")),
            all(vapply(fit$chains,function(ch) all(is.finite(ch$parameters)) &&
              all(is.finite(ch$beta)) && ncol(ch$beta)==2L,logical(1))))
  if (!global) stopifnot(fit$posterior_gp_included==1)
}
# No-design default and explicit mode choice must not alter the legacy draws.
small <- args; small$Z <- small$beta_prior_mean <- small$beta_prior_cov <- NULL
small$n_iter <- 30L; small$burnin <- 5L; small$n_chains <- 1L
f0 <- do.call(bmim_joint_domain_mcmc,small)
f1 <- do.call(bmim_joint_domain_mcmc,c(small,list(beta_mode="sampled")))
stopifnot(f0$beta_mode=="none",identical(f0,f1))
cat("SCALAR FULL-PARAMETER / WRAPPER / NO-DESIGN SMOKE: PASS\n")
