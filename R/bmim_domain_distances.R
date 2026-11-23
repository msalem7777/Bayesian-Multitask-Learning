# BMIM independent functional-domain geometry, priors, and proposals.
# Source bmim_joint_domains.R to load the complete joint-domain workflow.
# No GP or other package dependencies; all intervals use inclusive indices.

bmim_domain_prior <- function(T, max_length, min_length = 1L,
                              length_weights = NULL) {
  valid_int <- function(z) length(z) == 1L && is.numeric(z) &&
    is.finite(z) && z == as.integer(z)
  if (!valid_int(T) || !valid_int(max_length) ||
      !valid_int(min_length) || T < 1L || max_length < 1L ||
      min_length < 1L) {
    stop("T, max_length, and min_length must be positive integers.")
  }
  cap <- min(as.integer(T), as.integer(max_length))
  if (min_length > cap) {
    stop("min_length exceeds the available maximum domain length.")
  }
  lengths <- seq.int(as.integer(min_length), cap)
  if (is.null(length_weights)) length_weights <- rep(1, length(lengths))
  if (length(length_weights) != length(lengths) ||
      any(!is.finite(length_weights)) || any(length_weights <= 0)) {
    stop("length_weights must be strictly positive, one per allowed length.")
  }
  weights <- length_weights / sum(length_weights)
  list(T = as.integer(T), min_length = as.integer(min_length),
       max_length = cap, lengths = lengths, weights = weights)
}

bmim_domain_draw <- function(domain_prior) {
  k <- sample.int(length(domain_prior$lengths), 1L,
                  prob = domain_prior$weights)
  len <- domain_prior$lengths[k]
  start <- sample.int(domain_prior$T - len + 1L, 1L)
  c(start = start, end = start + len - 1L)
}

bmim_domain_log_prior <- function(A, domain_prior) {
  if (length(A) != 2L || any(!is.finite(A)) ||
      any(A != as.integer(A))) return(-Inf)
  start <- as.integer(A[1L])
  end <- as.integer(A[2L])
  len <- end - start + 1L
  k <- match(len, domain_prior$lengths)
  if (is.na(k) || start < 1L || end > domain_prior$T) return(-Inf)
  log(domain_prior$weights[k]) - log(domain_prior$T - len + 1L)
}

# Enumerate one-endpoint +/-1 and one-timepoint translation proposals.
# Translations preserve length, including when min_length=max_length.
# The neighborhood is symmetric as a graph, but degrees may differ at edges.
# The Hastings correction therefore uses log(degree_old / degree_new).
bmim_domain_neighbors <- function(A, domain_prior) {
  s <- as.integer(A[1L])
  e <- as.integer(A[2L])
  candidates <- rbind(c(s - 1L, e), c(s + 1L, e),
                      c(s, e - 1L), c(s, e + 1L),
                      c(s - 1L, e - 1L), c(s + 1L, e + 1L))
  keep <- apply(candidates, 1L, function(x) {
    is.finite(bmim_domain_log_prior(x, domain_prior))
  })
  out <- candidates[keep, , drop = FALSE]
  colnames(out) <- c("start", "end")
  out
}

bmim_domain_propose <- function(A, domain_prior,
                                prior_jump_probability = 0.50) {
  if (length(prior_jump_probability) != 1L ||
      !is.finite(prior_jump_probability) ||
      prior_jump_probability < 0 || prior_jump_probability > 1) {
    stop("prior_jump_probability must lie in [0,1].")
  }
  if (runif(1L) < prior_jump_probability) {
    # Independence proposal q(A*)=p(A*) exactly cancels domain prior.
    return(list(A = bmim_domain_draw(domain_prior),
                log_prior_and_proposal_ratio = 0, kind = "prior"))
  }
  old_neighbors <- bmim_domain_neighbors(A, domain_prior)
  if (nrow(old_neighbors) == 0L) {
    return(list(A = A, log_prior_and_proposal_ratio = 0,
                kind = "identity"))
  }
  new_A <- old_neighbors[sample.int(nrow(old_neighbors), 1L), ]
  new_neighbors <- bmim_domain_neighbors(new_A, domain_prior)
  if (nrow(new_neighbors) == 0L ||
      !any(new_neighbors[, 1L] == A[1L] & new_neighbors[, 2L] == A[2L])) {
    stop("Local domain proposal does not have a reverse move.")
  }
  correction <- bmim_domain_log_prior(new_A, domain_prior) -
    bmim_domain_log_prior(A, domain_prior) +
    log(nrow(old_neighbors)) - log(nrow(new_neighbors))
  list(A = new_A, log_prior_and_proposal_ratio = correction,
       kind = "local")
}

# Squared Euclidean between-subject distances averaged over one interval.
# Normalization by interval length matches the validated GP distance rule.
bmim_domain_distance <- function(X_s, loc_columns, A) {
  s <- as.integer(A[1L])
  e <- as.integer(A[2L])
  if (s < 1L || e > length(loc_columns) || s > e) {
    stop("The proposed interval lies outside the node's time grid.")
  }
  ids <- loc_columns[seq.int(s, e)]
  xx <- X_s[, ids, drop = FALSE] / sqrt(length(ids))
  norms <- rowSums(xx * xx)
  D <- outer(norms, norms, "+") - 2 * tcrossprod(xx)
  D <- (D + t(D)) / 2
  D[D < 0] <- 0
  diag(D) <- 0
  D
}

# Numerically stabilize a distance change made via incremental replacement.
bmim_joint_distance_replace <- function(total, old, new) {
  out <- total - old + new
  out <- (out + t(out)) / 2
  out[out < 0] <- 0
  diag(out) <- 0
  out
}

# The Ising and fixed VB empirical-Bayes node-selection prior.
# The empty configuration is forbidden, matching ISGM().
