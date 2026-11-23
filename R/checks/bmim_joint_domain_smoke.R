# ============================================================================
# bmim_joint_domain_smoke.R -- deterministic, independent reference checks.
# Run from R/:  source("checks/bmim_joint_domain_smoke.R")
# R dependencies: base R only.
# Tests: interval support, fixed-domain regression, proposal detailed balance,
# fixed-(d,g) exact enumeration incl. global b, and four-chain event estimates.
# ============================================================================
source_files <- Filter(Negate(is.null),
                       lapply(sys.frames(), function(frame) frame$ofile))
if (!length(source_files)) stop("Source this check from a real .R file.")
root <- dirname(dirname(normalizePath(tail(source_files, 1L)[[1L]])))
source(file.path(root, "bmim_gp_updates.R"))
source(file.path(root, "bmim_global_gp_test.R"))
source(file.path(root, "bmim_joint_domains.R"))
rm(source_files, root)

set.seed(20261122L)
n <- 7L
z <- seq(-1.4, 1.4, length.out = n)
X <- scale(cbind(z, z^2, sin(2*z),
                 cos(z), z + cos(2*z), sin(3*z)))
y <- as.numeric(scale(0.7*X[, 2L] + 0.4*X[, 5L] + 0.2*cos(4*z)))
loc <- list(1:3, 4:6)
W <- matrix(c(0, 0.3, 0.3, 0), 2L, 2L)
p_prop <- c(0.64, 0.38)
pi_b <- 0.50
prior <- bmim_gp_default_prior(bmim_gp_precompute_distances(X, loc))
fixed <- list(d = 0.6*prior$distance_scale, g = 0.65)

# Regression: if the only domain is the complete node, the new distance
# implementation is numerically identical to the previously validated one.
all_time <- bmim_domain_prior(T=3L, max_length=3L, min_length=3L)
stopifnot(all(bmim_domain_draw(all_time) == c(1L, 3L)))
for (p in 1:2) {
  old <- bmim_gp_precompute_distances(X, loc)[[p]]
  new <- bmim_domain_distance(X, loc[[p]], c(1L, 3L))
  stopifnot(isTRUE(all.equal(old, new, tolerance=1e-12)))
}
cat("FULL-DOMAIN GP DISTANCE REGRESSION: PASS\n")

# Five intervals per node when T=3 and max=2:
# 3 windows of length 1, 2 windows of length 2.
P <- bmim_domain_prior(T=3L, max_length=2L)
intervals <- rbind(c(1L,1L),c(2L,2L),c(3L,3L),c(1L,2L),c(2L,3L))
stopifnot(all(vapply(seq_len(nrow(intervals)), function(j)
  is.finite(bmim_domain_log_prior(intervals[j, ],P)), logical(1L))))
stopifnot(abs(sum(exp(apply(intervals,1L,bmim_domain_log_prior,
                            domain_prior=P))) - 1) < 1e-12)
stopifnot(all(vapply(seq_len(nrow(intervals)), function(j)
  length(bmim_domain_neighbors(intervals[j, ],P)[,1L]) >= 1L,
  logical(1L))))
cat("BOUNDED-DOMAIN PRIOR NORMALIZATION: PASS\n")

# Ensure the cached total D remains synchronized when the global indicator
# is OFF (b=0) but a node xi[p] is ON. The null likelihood ignores A, yet
# the next b-update still needs the GP-alternative covariance at the new A.
A_initial <- rbind(c(1L,1L),c(2L,2L))
xi_initial <- c(1L,0L)
D_initial <- lapply(1:2,function(p)
  bmim_domain_distance(X,loc[[p]],A_initial[p,]))
D_total_initial <- bmim_gp_selected_distance(D_initial,xi_initial)
update_under_null <- bmim_joint_update_domain(
  p=1L, xi=xi_initial, A=A_initial, D_by_loc=D_initial,
  D_total=D_total_initial, X_s=X, loc_slices=loc,
  domain_priors=list(P,P), y=y, b=0L, gp=list(d=1,g=0.5),
  gp_prior=prior, prior_jump_probability=0.5
)
D_rebuilt <- bmim_gp_selected_distance(update_under_null$D_by_loc,
                                        xi_initial)
stopifnot(isTRUE(all.equal(update_under_null$D_total,D_rebuilt,
                           tolerance=1e-10)))
cat("NULL-MODEL ACTIVE-DISTANCE CACHE CONSISTENCY: PASS\n")

# Direct independent full Gaussian-IG integrated likelihood, WITHOUT using
# the actual BMIM GP evidence routine. This is an oracle at fixed d,g.
ref_integrated <- function(y, covariance, shape, rate) {
  R <- chol(covariance)
  z <- forwardsolve(t(R), y)
  Q <- sum(z*z)
  logdet <- 2*sum(log(diag(R)))
  N <- length(y)
  lgamma(shape+N/2)-lgamma(shape)+shape*log(rate) -
    N/2*log(2*pi)-logdet/2 - (shape+N/2)*log(rate+Q/2)
}
ref_xiprior <- function(xi) {
  spins <- 2*xi-1
  sum(outer(spins,spins)*W)/2 +
    sum(xi*log(p_prop)+(1-xi)*log1p(-p_prop))
}
ref_domain_prior <- function(A) {
  L <- A[2L]-A[1L]+1L
  -log(2) - log(3-L+1L)  # uniform length in {1,2}, uniform start
}
ref_log_alt <- function(xi, A) {
  D <- matrix(0,n,n)
  for (p in which(xi==1L)) {
    cols <- loc[[p]][seq.int(A[p,1L],A[p,2L])]
    M <- X[,cols,drop=FALSE]/sqrt(length(cols))
    D <- D + outer(rowSums(M*M),rowSums(M*M),"+")-2*tcrossprod(M)
  }
  D <- (D+t(D))/2
  D[D<0] <- 0
  diag(D) <- 0
  C <- exp(-D/fixed$d)+diag(fixed$g,n)
  ref_integrated(y,C,prior$tau_shape,prior$tau_rate)
}
log_null <- ref_integrated(y,diag(n),prior$null_sigma_shape,
                           prior$null_sigma_rate)
stopifnot(isTRUE(all.equal(log_null,
  bmim_global_null_logevidence(y,prior),tolerance=1e-10)))

# All allowable states: three nonempty xi configurations and 5^2 domain
# combinations. The global indicator b adds the noise-only branch.
configs <- rbind(c(1L,0L),c(0L,1L),c(1L,1L))
states <- list()
logw <- numeric(0)
for (k in seq_len(nrow(configs))) {
  for (j1 in seq_len(nrow(intervals))) {
    for (j2 in seq_len(nrow(intervals))) {
      A <- rbind(intervals[j1, ],intervals[j2, ])
      xi <- configs[k, ]
      lp <- ref_xiprior(xi) +
        ref_domain_prior(A[1L, ]) + ref_domain_prior(A[2L, ])
      log_alt <- ref_log_alt(xi,A)
      for (b in 0:1) {
        states[[length(states)+1L]] <- list(b=b, xi=xi, A=A)
        log_branch <- if (b == 1L) {
          log(pi_b) + log_alt
        } else {
          log1p(-pi_b) + log_null
        }
        logw <- c(logw, lp + log_branch)
      }
    }
  }
}
weights <- exp(logw-max(logw))
weights <- weights/sum(weights)
stopifnot(length(states)==150L,abs(sum(weights)-1)<1e-12)
exact_b <- sum(weights[vapply(states,function(s) s$b==1L,logical(1L))])
exact_nodes <- vapply(1:2,function(p) sum(weights[vapply(states,function(s)
  s$xi[p]==1L,logical(1L))]),numeric(1L))
exact_joint <- matrix(0,2L,3L)
for (p in 1:2) for (t in 1:3) {
  exact_joint[p,t] <- sum(weights[vapply(states,function(s)
    s$b==1L && s$xi[p]==1L &&
    s$A[p,1L]<=t && s$A[p,2L]>=t,logical(1L))])
}
cat("INDEPENDENT 150-STATE POSTERIOR REFERENCE: PASS\n")

# Verify the node+domain joint proposal reversibility algebra, including
# the p(A_new)/p(A_old) terms which cancel the independence proposal ratio.
for (k in 1:3) {
  xi <- configs[k,]
  for (p in 1:2) {
    xi_new <- xi
    xi_new[p] <- 1L-xi_new[p]
    if (sum(xi_new)==0L) next
    for (j1 in 1:5) for (j2 in 1:5) {
      old_A <- rbind(intervals[j1,],intervals[j2,])
      for (jstar in 1:5) {
        new_A <- old_A
        new_A[p,] <- intervals[jstar,]
        lp_old <- ref_xiprior(xi) + ref_domain_prior(old_A[p,]) +
          ref_log_alt(xi,old_A)
        lp_new <- ref_xiprior(xi_new) + ref_domain_prior(new_A[p,]) +
          ref_log_alt(xi_new,new_A)
        q_forward <- ref_domain_prior(new_A[p,])
        q_reverse <- ref_domain_prior(old_A[p,])
        lr <- ref_xiprior(xi_new) - ref_xiprior(xi) +
          ref_log_alt(xi_new,new_A) - ref_log_alt(xi,old_A)
        flow_forward <- lp_old+q_forward+min(0,lr)
        flow_reverse <- lp_new+q_reverse+min(0,-lr)
        stopifnot(isTRUE(all.equal(flow_forward,flow_reverse,
                                   tolerance=1e-9)))
      }
    }
  }
}
cat("JOINT NODE/DOMAIN PROPOSAL DETAILED BALANCE: PASS\n")

# Detailed balance for local endpoint moves, including boundary proposal
# probabilities which are non-uniform because neighbor counts differ.
for (i in 1:5) {
  Ai <- intervals[i,]
  neighbors <- bmim_domain_neighbors(Ai,P)
  for (j in seq_len(nrow(neighbors))) {
    Aj <- neighbors[j,]
    reverse_neighbors <- bmim_domain_neighbors(Aj,P)
    li <- bmim_domain_log_prior(Ai,P)
    lj <- bmim_domain_log_prior(Aj,P)
    lr <- lj-li+log(nrow(neighbors))-log(nrow(reverse_neighbors))
    flow_forward <- li-log(nrow(neighbors))+min(0,lr)
    flow_reverse <- lj-log(nrow(reverse_neighbors))+min(0,-lr)
    stopifnot(isTRUE(all.equal(flow_forward,flow_reverse,
                               tolerance=1e-12)))
  }
}
cat("LOCAL DOMAIN PROPOSAL DETAILED BALANCE: PASS\n")

# A discrete MCMC with fixed d,g must reproduce the independently enumerated
# posterior including the global Bernoulli indicator. Four different chains.
fit <- bmim_joint_domain_mcmc(
  X_s=X,y_tilde=y,loc_slices=loc,W_mat=W,p_prop=p_prop,
  max_domain_length=2L,n_iter=3500L,burnin=700L,n_chains=4L,
  min_domain_length=1L,pi_b=pi_b,
  n_domain_updates_per_sweep=1L, domain_prior_jump_prob=0.5,
  gp_prior=prior,fixed_gp=fixed,seed=20261122L
)
stopifnot(abs(fit$posterior_gp_included-exact_b)<0.045)
stopifnot(max(abs(fit$posterior_node_included-exact_nodes))<0.055)
stopifnot(max(abs(fit$posterior_node_time_joint_gp-exact_joint))<0.055)
for (chain in fit$chains) {
  stopifnot(all(chain$end>=chain$start),
            all(chain$end-chain$start+1L<=2L),
            all(rowSums(chain$xi)>=1L),
            all(is.finite(chain$parameters)),
            all(chain$parameters[,"d"]==fixed$d),
            all(chain$parameters[,"g"]==fixed$g))
}
cat("FIXED-D/G JOINT DOMAIN+GLOBAL BERNOULLI MCMC: PASS\n")
cat(sprintf("P(b=1): exact %.5f / sampled %.5f\n",
            exact_b,fit$posterior_gp_included))
print(data.frame(node=1:2,exact=round(exact_nodes,4),
                 sampled=round(fit$posterior_node_included,4)))
cat("BMIM JOINT DOMAIN SMOKE: PASS\n")
