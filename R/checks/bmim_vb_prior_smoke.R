# Dedicated VB tests. Base-R algebra plus actual ClusterR/infotheo integration.
# Run from R/: Rscript --vanilla -e "source('checks/bmim_vb_prior_smoke.R')"
files <- Filter(Negate(is.null),lapply(sys.frames(),function(f) f$ofile))
root <- dirname(dirname(normalizePath(tail(files,1L)[[1L]])))
source(file.path(root,"bmim_gp_updates.R"))
source(file.path(root,"bmim_global_gp_test.R"))
source(file.path(root,"bmim_joint_domains.R"))
source(file.path(root,"ISGM_core.R"))
source(file.path(root,"bmim_joint_simulation.R"))

# Compare the stable log ratio to the algebraic Gaussian density identity,
# including values at which the old direct density ratio underflowed.
u <- c(.01,.05,.2,.4)
expected <- -log(2)+.5*log(.03/.04)+u^2/(2*.03)-(u-.2)^2/(2*.04)
stopifnot(max(abs(bmim_vb_logodds(u,.2,.04,.03)-expected))<1e-12,
          all(is.finite(bmim_vb_logodds(c(.1,.5),.2,1e-8,1e-8))))
cat("VB STABLE WORKING LOG DENSITIES: PASS\n")

# Independently enumerate the variational expectation of the unnormalized
# auxiliary log distribution and entropy, for all 2^P binary configurations.
delta <- c(-.6,.4,-.2,.7)
W <- matrix(c(0,.3,-.2,.1,.3,0,.4,.2,-.2,.4,0,.15,.1,.2,.15,0),4)
ref_objective <- function(m) {
  x <- as.matrix(expand.grid(rep(list(0:1),length(m))))
  pp <- (1+m)/2
  q <- apply(x,1,function(z) prod(ifelse(z==1,pp,1-pp)))
  s <- 2*x-1
  logp <- as.numeric(x%*%delta)+rowSums((s%*%W)*s)/2
  keep <- q>0
  sum(q[keep]*(logp[keep]-log(q[keep])))
}
set.seed(991)
one <- bmim_vb_mean_field(delta,W,iters=1L)
set.seed(991)
ref <- tanh(delta/2)
for (j in sample.int(4L,replace=FALSE)) {
  # Conditional odds computed independently from the two coordinate states.
  plus <- ref; minus <- ref; plus[j]<-1; minus[j]<- -1
  logratio <- ref_objective(plus)-ref_objective(minus)
  ref[j] <- 2*plogis(logratio)-1
}
stopifnot(max(abs(ref-one$means))<1e-12, all(one$updates==1L),
          abs(ref_objective(one$means)-tail(one$objective,1))<1e-12)
set.seed(992)
fit <- bmim_vb_mean_field(delta,W,iters=15L)
stopifnot(all(fit$updates==15L),all(diff(fit$objective)> -1e-10),
          abs(ref_objective(fit$means)-tail(fit$objective,1))<1e-12,
          max(abs(bmim_vb_mean_field(delta,0*W)$weights-plogis(delta)))<1e-12)
cat("VB FULL COORDINATE SWEEPS / EXACT ELBO REFERENCE: PASS\n")

# Scalar confounding: X1 shares Z with y, X2 retains actual residual signal.
# Construct orthogonal vectors so the expected partial correlations are exact.
n <- 180L; t <- seq_len(n)
z <- cos(2*pi*t/n); a <- sin(2*pi*t/n); e <- sin(4*pi*t/n)
y <- 5*z+a+e; X <- cbind(5*z+e*0+sin(6*pi*t/n),a,z)
loc <- as.list(seq_len(ncol(X))); Z <- cbind(z,z,1)
raw <- bmim_vb_associations(X,y,loc,vb_summary="sq.corr")
adj <- bmim_vb_associations(X,y,loc,Z,"partial.sq.corr")
stopifnot(raw$node[1]>.8,adj$node[1]<1e-12,
          abs(adj$node[2]-.5)<1e-12,adj$node[3]==0,
          identical(unname(adj$zero_columns),3L),
          bmim_vb_method("auto")=="sq.corr",
          bmim_vb_method("auto",Z)=="partial.sq.corr",
          bmim_vb_method()=="sq.corr",
          bmim_vb_method(Z=Z)=="partial.sq.corr",
          bmim_vb_method("sq.corr",Z)=="sq.corr",
          bmim_vb_method("mi",Z)=="mi",
          bmim_vb_method("residual.mi",Z)=="residual.mi")
# Default association calculation must match the explicit correlation choice.
stopifnot(identical(bmim_vb_associations(X,y,loc),raw),
          identical(bmim_vb_associations(X,y,loc,Z),adj))
stopifnot(inherits(try(bmim_vb_method("conditional.mi",Z),silent=TRUE),"try-error"),
          inherits(try(bmim_vb_method("residual.mi"),silent=TRUE),"try-error"))
neutral <- bmim_vb_fit(rep(.1,4),W)
stopifnot(all(neutral$weights==.5),neutral$status=="identical_summaries_neutral_evidence")
cat("VB SCALAR ADJUSTMENT / DEFAULTS / DEGENERATE SUMMARIES: PASS\n")

# Package-backed checks are mandatory, not silently skipped.
for (pkg in c("ClusterR","infotheo"))
  if (!requireNamespace(pkg,quietly=TRUE)) stop("VB smoke requires package: ",pkg)
# Independent historical discretized empirical MI calculation (same explicit bins).
# Quadratic dependence should be much stronger than squared Pearson correlation.
x <- seq(-1,1,length.out=n); yy <- x^2+.01*cos(2*pi*t/n)
XX <- cbind(x, sin(17*pi*t/n)); ll <- list(1L,2L)
mi <- bmim_vb_associations(XX,yy,ll,vb_summary="mi")
co <- bmim_vb_associations(XX,yy,ll,vb_summary="sq.corr")
ref_mi <- infotheo::mutinformation(infotheo::discretize(yy-mean(yy),nbins=mi$bins),
                                   infotheo::discretize(x-mean(x),nbins=mi$bins),method="emp")
stopifnot(abs(mi$node[1]-ref_mi)<1e-12,mi$node[1]>.1,co$node[1]<.01)
res <- bmim_vb_associations(X,y,loc,Z,"residual.mi")
stopifnot(res$node[3]==0,all(is.finite(res$node)),res$method=="residual.mi")
# Use clearly separated nondegenerate data for the actual GMM path, and
# ensure input containers cannot change which coordinates are updated.
u <- c(.01,.02,.035,.05,.065,.08,.24,.27,.30,.34,.38,.43)
WW <- matrix(.035,12,12); diag(WW)<-0
set.seed(993); a1 <- bmim_vb_fit(matrix(u,ncol=1),WW)
set.seed(993); a2 <- bmim_vb_fit(data.frame(u),WW)
stopifnot(isTRUE(all.equal(a1$weights,a2$weights,tolerance=1e-12)),
          all(a1$updates==15L),all(diff(a1$objective)> -1e-9))
set.seed(993); old_interface <- VBprop(data.frame(u),X=NULL,corrmat=WW,num_locs=12,iters=15)
stopifnot(isTRUE(all.equal(old_interface,a1$weights,tolerance=1e-12)))
cat("VB ACTUAL GMM / MI / DATA-FRAME REGRESSION: PASS\n")

# Feature-wrapper integration with all four measures. Explicit weights must
# bypass VB, and automatic weights are generated exactly once before MCMC.
set.seed(994)
n <- 90L; z <- rnorm(n); y <- .6*z+rnorm(n)
X <- sapply(seq(.02,.9,length.out=12),function(r) r*y+sqrt(1-r*r)*rnorm(n))
loc <- as.list(seq_len(ncol(X))); W <- matrix(.025,12,12); diag(W)<-0
for (method in c("sq.corr","partial.sq.corr","mi","residual.mi")) {
  Z <- if (method %in% c("partial.sq.corr","residual.mi")) matrix(z,ncol=1) else NULL
  args <- list(path_features=X,W_mat=W,y=y,grps_vec=loc,max_domain_length=1L,
               Z=Z,vb_summary=method,vb_seed=995L,n_iter=8L,burnin=2L,n_chains=1L,
               fixed_gp=list(d=1,g=.3),pi_b=.5,seed=996L)
  ans <- do.call(bmim_joint_domains_from_features,args)
  stopifnot(ans$vb$association$method==method,all(ans$vb$updates==15L),
            identical(ans$p_prop,ans$vb$weights),is.finite(ans$posterior_gp_included))
  args$p_prop <- rep(.5,12)
  explicit <- do.call(bmim_joint_domains_from_features,args)
  stopifnot(is.null(explicit$vb),all(explicit$p_prop==.5))
}
# Verify the default wrapper dispatch, with and without scalar covariates.
for (has_z in c(FALSE,TRUE)) {
  Z_default <- if (has_z) matrix(z,ncol=1) else NULL
  args <- list(path_features=X,W_mat=W,y=y,grps_vec=loc,max_domain_length=1L,
               Z=Z_default,vb_seed=995L,n_iter=8L,burnin=2L,n_chains=1L,
               fixed_gp=list(d=1,g=.3),pi_b=.5,seed=996L)
  automatic <- do.call(bmim_joint_domains_from_features,args)
  args$vb_summary <- if (has_z) "partial.sq.corr" else "sq.corr"
  explicit <- do.call(bmim_joint_domains_from_features,args)
  stopifnot(automatic$vb$association$method==args$vb_summary,
            identical(automatic$p_prop,explicit$p_prop),
            identical(automatic$posterior_gp_included,explicit$posterior_gp_included))
}
# Simulation adapter uses the same constructor and keeps flat mode unchanged.
set.seed(997); old_rng <- .Random.seed
v <- bmim_simulation_selection_prior(X,y,loc,W,"vb",vb_summary="mi",details=TRUE,vb_seed=998L)
stopifnot(identical(old_rng,.Random.seed),v$association$method=="mi",all(v$updates==15L),
          all(bmim_simulation_selection_prior(X,y,loc,W,"flat")==.5))
default_v <- bmim_simulation_selection_prior(X,y,loc,W,"vb",details=TRUE,vb_seed=998L)
explicit_v <- bmim_simulation_selection_prior(X,y,loc,W,"vb",vb_summary="sq.corr",details=TRUE,vb_seed=998L)
stopifnot(default_v$association$method=="sq.corr",identical(default_v$weights,explicit_v$weights))
cat("VB CORRELATION DEFAULT DISPATCH: PASS\n")
cat("VB FEATURE WRAPPER / SIMULATION ADAPTER / RNG: PASS\n")
cat("BMIM VB PRIOR SMOKE: PASS\n")
