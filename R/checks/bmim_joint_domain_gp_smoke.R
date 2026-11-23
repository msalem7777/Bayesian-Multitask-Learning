# ============================================================================
# bmim_joint_domain_gp_smoke.R
# Short full-parameter integration smoke: d,g sampled, node/intervals/b moved.
# This is NOT a convergence or statistical-calibration study.
# Run from R/:  source("checks/bmim_joint_domain_gp_smoke.R")
# ============================================================================
files <- Filter(Negate(is.null),
                lapply(sys.frames(), function(frame) frame$ofile))
if (!length(files)) stop("Source this test as a file.")
root <- dirname(dirname(normalizePath(tail(files,1L)[[1L]])))
source(file.path(root,"bmim_gp_updates.R"))
source(file.path(root,"bmim_global_gp_test.R"))
source(file.path(root,"bmim_joint_domains.R"))
rm(files,root)

set.seed(20261122L)
z <- seq(-1.2,1.2,length.out=8L)
X <- scale(cbind(z,z^2,sin(z),cos(z),
                 cos(2*z),z^3,z+sin(z),sin(3*z)))
y <- as.numeric(scale(0.4*X[,2L]+0.7*X[,5L]+0.1*sin(4*z)))
loc <- list(1:4,5:8)
W <- matrix(c(0,0.25,0.25,0),2L,2L)
p_prop <- c(0.65,0.42)

# Test selection-only mode first: no global b, with explicit max length.
node_only <- bmim_joint_domain_mcmc(
  X,y,loc,W,p_prop,max_domain_length=2L,
  n_iter=80L,burnin=20L,n_chains=2L,pi_b=NULL,seed=20261122L
)
stopifnot(node_only$posterior_gp_included==1,
          all(is.finite(node_only$posterior_node_time_joint_gp)),
          all(node_only$posterior_node_time_joint_gp>=0),
          all(node_only$posterior_node_time_joint_gp<=1))

# Test actual GP hyperparameter updates + b switching together.
global <- bmim_joint_domain_mcmc(
  X,y,loc,W,p_prop,max_domain_length=3L,
  n_iter=120L,burnin=30L,n_chains=2L,pi_b=0.45,
  n_domain_updates_per_sweep=1L,seed=20261123L
)
for (chain in global$chains) {
  stopifnot(all(chain$end>=chain$start),
            all(chain$end-chain$start+1L<=3L),
            all(rowSums(chain$xi)>0),
            all(chain$parameters[,"d"]>0),
            all(chain$parameters[,"g"]>0),
            all(chain$parameters[,"tau2"]>0),
            all(chain$parameters[,"sigma2"]>0),
            all(is.finite(chain$parameters)))
}
stopifnot(is.finite(global$posterior_gp_included),
          global$posterior_gp_included>=0,
          global$posterior_gp_included<=1,
          all(global$posterior_node_time_joint_gp>=0,na.rm=TRUE),
          all(global$posterior_node_time_joint_gp<=1,na.rm=TRUE))
cat("BMIM JOINT DOMAIN FULL-PARAMETER INTEGRATION SMOKE: PASS\n")
