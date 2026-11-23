# ============================================================================
# Eight-state global GP indicator validation (base R only)
#
# Enumerates one exact noise-only state and seven nonempty GP configurations.
# GP parameters are held fixed to make the eight-state posterior analytically
# computable without numerical integration. Four independent Gibbs chains
# are compared with the exact reference frequencies. This complements the
# full-parameter sampler smoke test; it does not establish calibration of the
# response-adaptive empirical-Bayes model.
# ============================================================================

.test_sources <- Filter(Negate(is.null),
                        lapply(sys.frames(), function(frame) frame$ofile))
if (!length(.test_sources)) stop('Source this test from its .R file.')
.test_file <- normalizePath(tail(.test_sources, 1L)[[1L]])
.test_root <- dirname(dirname(.test_file))
source(file.path(.test_root, 'bmim_gp_updates.R'))
source(file.path(.test_root, 'bmim_global_gp_test.R'))
source(file.path(.test_root, 'bmim_global_gp_mcmc.R'))
rm(.test_sources, .test_file, .test_root)

set.seed(20261122L)
n <- 9L
z <- seq(-1.1, 1.1, length.out=n)
X <- scale(cbind(z, sin(z), cos(2*z), z^2, sin(3*z), z^3))
y <- as.numeric(scale(0.7*X[,1] - 0.5*X[,3] + 0.2*sin(4*z)))
D <- bmim_gp_precompute_distances(X, list(1:2, 3:4, 5:6))
prior <- bmim_gp_default_prior(D)
configs <- as.matrix(expand.grid(rep(list(0:1), 3L)))
configs <- configs[rowSums(configs) > 0L, , drop=FALSE]
storage.mode(configs) <- 'double'
stopifnot(nrow(configs) == 7L)
W <- matrix(c(0, .25, .1, .25, 0, .2, .1, .2, 0), nrow=3L)
r <- c(.7, .4, .6)
pi_b <- .4
d <- prior$distance_scale
g <- .7
tau2 <- 1.2

# Independent eight-state reference with both variance scales integrated.
# Under b=0: sigma0^2~IG(a0,r0), covariance sigma0^2 I.
# Under b=1: tau2~IG(a1,r1), covariance tau2 (K_xi(d)+gI).
independent_integrated_gaussian <- function(y, C, shape, rate) {
  R <- chol(C)
  v <- forwardsolve(t(R), y)
  Q <- sum(v*v)
  ld <- 2*sum(log(diag(R)))
  n <- length(y)
  lgamma(shape+n/2) - lgamma(shape) + shape*log(rate) -
    n/2*log(2*pi) - ld/2 - (shape+n/2)*log(rate+Q/2)
}
log_prior <- vapply(seq_len(nrow(configs)), function(k) {
  xi <- configs[k, ]
  spin <- 2*xi - 1
  as.numeric(crossprod(spin, W %*% spin))/2 +
    sum(xi*log(r)+(1-xi)*log1p(-r))
}, numeric(1))
log_prior <- log_prior - max(log_prior)
model_prior <- exp(log_prior) / sum(exp(log_prior))
log_null <- independent_integrated_gaussian(
  y, diag(n), prior$null_sigma_shape, prior$null_sigma_rate
)
log_alt <- vapply(seq_len(nrow(configs)), function(k) {
  Dk <- bmim_gp_selected_distance(D, configs[k, ])
  Ck <- exp(-Dk/d)+g*diag(n)
  independent_integrated_gaussian(
    y, Ck, prior$tau_shape, prior$tau_rate
  )
}, numeric(1))
# Cross-check against the production marginal likelihood routines.
stopifnot(isTRUE(all.equal(
  log_null, bmim_global_null_logevidence(y, prior), tolerance=1e-10)))
for (k in seq_len(nrow(configs))) {
  Dk <- bmim_gp_selected_distance(D, configs[k, ])
  helper <- bmim_global_log_integrated_tau(
    y, bmim_gp_covariance(Dk, d, g), prior
  )
  stopifnot(isTRUE(all.equal(as.numeric(helper), as.numeric(log_alt[k]),
                             tolerance=1e-10)))
}

# Eight mutually exclusive states: null=0, GP with each configuration=1..7.
log_joint <- c(log1p(-pi_b)+log_null,
               log(pi_b)+log(model_prior)+log_alt)
weights <- exp(log_joint-max(log_joint))
exact <- weights/sum(weights)
stopifnot(length(exact)==8L, abs(sum(exact)-1)<1e-12)

# Direct Gibbs transitions at FIXED d,g, with both variance scales integrated. They constitute an independent
# test of the joint indicator/configuration target, not a replacement for the
# full GP-parameter MCMC test in bmim_global_gp_smoke.R.
run_chain <- function(seed, n_iter=6000L, burnin=1000L) {
  set.seed(seed)
  b <- as.integer(seed %% 2L)
  idx <- sample.int(7L, 1L, prob=model_prior)
  states <- integer(n_iter-burnin)
  for (it in seq_len(n_iter)) {
    # Draw the configuration from its exact full conditional.
    if (b == 0L) {
      idx <- sample.int(7L, 1L, prob=model_prior)
    } else {
      q <- exp(log(model_prior)+log_alt-
               max(log(model_prior)+log_alt))
      idx <- sample.int(7L, 1L, prob=q)
    }
    # Draw the global indicator from its exact full conditional.
    prob_b <- plogis(qlogis(pi_b)+log_alt[idx]-log_null)
    b <- rbinom(1L, 1L, prob_b)
    if (it > burnin) states[it-burnin] <- if (b==0L) 1L else idx+1L
  }
  tabulate(states, nbins=8L)/length(states)
}
chain_prob <- do.call(rbind, lapply(1:4, function(j) run_chain(1000L+j)))
empirical <- colMeans(chain_prob)
error <- abs(empirical-exact)
cat('EIGHT-STATE REFERENCE AND FOUR-CHAIN FREQUENCIES:\n')
print(data.frame(state=c('null',apply(configs,1L,paste0,collapse='')),
                 exact=round(exact,5),
                 empirical=round(empirical,5),
                 absolute_error=round(error,5)), row.names=FALSE)
cat(sprintf('Global GP probability: exact %.5f, sampled %.5f\n',
            1-exact[1L],1-empirical[1L]))
stopifnot(max(error)<0.04,
          max(abs(sweep(chain_prob, 2L, exact, '-')))<0.07,
          abs((1-exact[1L])-(1-empirical[1L]))<0.04)
cat('EIGHT-STATE GLOBAL-INDICATOR TOY: PASS\n')
