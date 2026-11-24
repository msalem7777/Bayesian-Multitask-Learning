# Standalone empirical Ising fields: association summaries, working densities,
# and coordinate-ascent mean field. No GP likelihood is modified here.
bmim_vb_method <- function(vb_summary = "auto", Z = NULL) {
  method <- match.arg(vb_summary, c("auto", "sq.corr", "partial.sq.corr",
                                     "mi", "residual.mi"))
  if (method == "auto") method <- if (is.null(Z)) "mi" else "residual.mi"
  if (method %in% c("partial.sq.corr", "residual.mi") && is.null(Z))
    stop("Scalar-adjusted VB summaries require Z (use the unadjusted variant otherwise).")
  method
}

# X contains complete trajectories; loc_slices partitions its columns.
# Z is the scalar design (an intercept is included here even if not supplied).
# QR removes the linear span without explicitly forming the projection matrix.
bmim_vb_associations <- function(X, y, loc_slices, Z = NULL,
                                 vb_summary = "auto") {
  X <- as.matrix(X)
  if (!is.numeric(X) || !is.numeric(y) || nrow(X) < 3L ||
      length(y) != nrow(X) || any(!is.finite(X)) || any(!is.finite(y)) ||
      !length(loc_slices) || any(lengths(loc_slices) < 1L) ||
      !identical(sort(as.integer(unlist(loc_slices))), seq_len(ncol(X))))
    stop("Invalid complete functional data or node slices for VB.")
  y <- as.numeric(y)
  method <- bmim_vb_method(vb_summary, Z)
  # Center first for numerical stability; test residual variation relative to
  # original centered variation, so the rule is invariant to measurement units.
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  xnorm <- sqrt(colSums(Xc^2)); ynorm <- sqrt(sum(yc^2))
  if (method %in% c("partial.sq.corr", "residual.mi")) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) || nrow(Z) != nrow(X) || any(!is.finite(Z)))
      stop("Z must be a finite numeric scalar design with one row per subject.")
    # Adding an intercept twice is unnecessary; center the supplied columns,
    # drop constants, scale nonconstant columns and then prepend an intercept.
    Zc <- scale(Z, center = TRUE, scale = FALSE)
    zn <- sqrt(colSums(Zc^2)); keep <- zn > 0
    if (any(keep)) Zc <- sweep(Zc[, keep, drop=FALSE], 2L, zn[keep], "/")
    else Zc <- matrix(numeric(0), nrow(X), 0L)
    qrZ <- qr(cbind(1, Zc))
    if (nrow(X) - qrZ$rank < 3L)
      stop("Fewer than three residual degrees of freedom for scalar-adjusted VB.")
    Xc <- qr.resid(qrZ, Xc)
    yc <- as.numeric(qr.resid(qrZ, yc))
  }
  # A constant or numerically annihilated variable carries no residual
  # association evidence. Define its score as zero (correlation is undefined).
  tol <- sqrt(.Machine$double.eps)
  good <- xnorm > 0 & sqrt(colSums(Xc^2)) > tol*xnorm
  good_y <- ynorm > 0 && sqrt(sum(yc^2)) > tol*ynorm
  scores <- numeric(ncol(X))
  use_mi <- method %in% c("mi", "residual.mi")
  if (use_mi && !requireNamespace("infotheo", quietly = TRUE))
    stop("MI VB summaries require infotheo; install it or explicitly select a correlation summary.")
  if (good_y && any(good)) {
    if (use_mi) {
      # Restore the historical empirical, equal-frequency discretized MI.
      # Explicit bins match infotheo's historical default N^(1/3), rounded down.
      bins <- max(2L, as.integer(floor(nrow(X)^(1/3))))
      yd <- infotheo::discretize(yc, disc = "equalfreq", nbins = bins)
      for (j in which(good)) {
        xd <- infotheo::discretize(Xc[,j], disc = "equalfreq", nbins = bins)
        scores[j] <- max(0, infotheo::mutinformation(yd, xd, method = "emp"))
      }
    } else {
      scores[good] <- as.numeric(cor(yc, Xc[,good,drop=FALSE])^2)
    }
  }
  if (any(!is.finite(scores))) stop("Nonfinite VB association summaries.")
  list(method = method, column = scores,
       node = vapply(loc_slices, function(idx) mean(scores[idx]), numeric(1L)),
       zero_columns = which(!good), zero_response = !good_y,
       bins = if (use_mi) max(2L, as.integer(floor(nrow(X)^(1/3)))) else NULL)
}

# Gaussian inclusion and half-normal exclusion; calculate on log scale.
bmim_vb_logodds <- function(u, mu_plus, var_plus, var_minus) {
  if (any(!is.finite(c(u, mu_plus, var_plus, var_minus))) ||
      any(u < 0) || length(mu_plus) != 1L || length(var_plus) != 1L ||
      length(var_minus) != 1L || var_plus <= 0 || var_minus <= 0)
    stop("Invalid association values or working mixture parameters.")
  ans <- dnorm(u, mu_plus, sqrt(var_plus), log=TRUE) -
    dnorm(u, 0, sqrt(var_minus), log=TRUE) - log(2)
  if (any(!is.finite(ans))) stop("Nonfinite VB log evidence ratios.")
  ans
}

# ELBO up to an additive constant independent of Q: use Delta, not a
# second logistic transform of the current means. The Ising normalizer and
# exclusion log densities are constants for this auxiliary optimization.
bmim_vb_objective <- function(m, delta, W) {
  p <- (1+m)/2
  ent <- function(x) { ans <- numeric(length(x)); ix <- x > 0
    ans[ix] <- -x[ix]*log(x[ix]); ans }
  sum(p*delta) + 0.5*sum(m*as.numeric(W%*%m)) + sum(ent(p)+ent(1-p))
}

bmim_vb_mean_field <- function(delta, W, iters = 15L, lambda = 1,
                               clip = 1e-8) {
  delta <- as.numeric(delta); P <- length(delta); W <- as.matrix(W)
  if (!P || any(!is.finite(delta)) || !is.numeric(W) ||
      !identical(dim(W), c(P,P)) || any(!is.finite(W)) ||
      max(abs(W-t(W))) > 1e-9 || any(abs(diag(W)) > 1e-9) ||
      length(iters)!=1L || !is.finite(iters) || iters < 1L || iters != floor(iters) ||
      length(lambda)!=1L || !is.finite(lambda) || lambda <= 0 || lambda > 1 ||
      length(clip)!=1L || !is.finite(clip) || clip <= 0 || clip >= .5)
    stop("Invalid mean-field inputs; W must be symmetric with zero diagonal.")
  W <- (W+t(W))/2; diag(W) <- 0
  m <- tanh(delta/2)
  objective <- numeric(iters+1L); objective[1L] <- bmim_vb_objective(m,delta,W)
  updates <- integer(P); max_change <- numeric(iters)
  for (i in seq_len(iters)) {
    old <- m
    # P, not length(data.frame), is the number of nodes.
    for (j in sample.int(P, replace=FALSE)) {
      m[j] <- (1-lambda)*m[j] + lambda*tanh(sum(W[j,]*m)+delta[j]/2)
      updates[j] <- updates[j]+1L
    }
    objective[i+1L] <- bmim_vb_objective(m,delta,W)
    max_change[i] <- max(abs(m-old))
  }
  list(weights=pmin(pmax((1+m)/2,clip),1-clip), means=m,
       logodds=delta, objective=objective, updates=updates, max_change=max_change)
}

# Fit the historical Gaussian mixture, retaining the higher component mean
# and both variances; fitted mixture proportions are NOT used in Delta.
bmim_vb_fit <- function(u, W, iters=15L, lambda=1, clip=1e-8) {
  if (!(is.numeric(u) || is.data.frame(u)) || (is.data.frame(u) && ncol(u)!=1L) ||
      (is.matrix(u) && ncol(u)!=1L)) stop("u must be one association value per node.")
  u <- as.numeric(as.matrix(u))
  if (!length(u) || any(!is.finite(u)) || any(u<0)) stop("Invalid VB association vector.")
  # No empirical separation can be estimated in the degenerate case. Use
  # neutral external fields, not an arbitrary two-component mixture. Record it.
  if (length(unique(u)) < 2L) {
    out <- bmim_vb_mean_field(rep(0,length(u)),W,iters,lambda,clip)
    out$mixture <- NULL; out$status <- "identical_summaries_neutral_evidence"
    return(out)
  }
  if (!requireNamespace("ClusterR", quietly=TRUE)) stop("VB mixture requires ClusterR.")
  mix <- ClusterR::GMM(matrix(u,ncol=1L), gaussian_comps=2,
    dist_mode="eucl_dist", seed_mode="random_subset", km_iter=10, em_iter=50)
  means <- as.numeric(mix$centroids); vars <- as.numeric(mix$covariance_matrices)
  if (length(means)!=2L || length(vars)!=2L || any(!is.finite(c(means,vars))) || any(vars<=0))
    stop("Invalid fitted VB mixture; inspect association summaries rather than silently substituting a prior.")
  ord <- order(means); low <- ord[1L]; high <- ord[2L]
  delta <- bmim_vb_logodds(u,means[high],vars[high],vars[low])
  out <- bmim_vb_mean_field(delta,W,iters,lambda,clip)
  out$mixture <- list(mu_minus=means[low],mu_plus=means[high],
                      var_minus=vars[low],var_plus=vars[high])
  out$status <- "fitted"
  out
}

# Public constructor. Optionally isolate its RNG from the caller.
bmim_vb_prior <- function(X,y,loc_slices,W,Z=NULL,vb_summary="auto",
                          iters=15L,seed=NULL) {
  if (!is.null(seed)) {
    if (length(seed)!=1L || !is.finite(seed) || seed<1 || seed>.Machine$integer.max || seed!=floor(seed))
      stop("vb_seed must be a positive integer.")
    had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    if (had) saved <- get(".Random.seed",envir=.GlobalEnv)
    on.exit(if (had) assign(".Random.seed",saved,envir=.GlobalEnv) else
      if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
    set.seed(as.integer(seed))
  }
  assoc <- bmim_vb_associations(X,y,loc_slices,Z,vb_summary)
  out <- bmim_vb_fit(assoc$node,W,iters=iters)
  out$association <- assoc; out$version <- "vb2"
  out
}
