# BMIM reversible active-interval and joint node/interval MH kernels.
# Requires bmim_domain_distances.R and validated bmim_gp_updates.R.
# Source bmim_joint_domains.R to load the complete workflow.

bmim_joint_log_xi_prior <- function(xi, W_mat, p_prop) {
  if (length(xi) != length(p_prop) || any(!xi %in% c(0, 1)) ||
      sum(xi) == 0L) return(-Inf)
  spins <- 2 * xi - 1
  0.5 * as.numeric(crossprod(spins, W_mat %*% spins)) +
    sum(dbinom(xi, size = 1L, prob = p_prop, log = TRUE))
}

# This likelihood differs from the full normalizing-constant version only
# by terms independent of xi/A/d/g; use it for within-model MH proposals.
bmim_joint_log_alt_kernel <- function(y, D_total, gp, gp_prior) {
  bmim_gp_log_marginal_likelihood(y, D_total, gp$d, gp$g, gp_prior)
}

# Update one active domain. A prior-independent proposal or a symmetric-graph
# local proposal leaves the exact target invariant after the MH correction.
bmim_joint_update_domain <- function(p, xi, A, D_by_loc, D_total,
                                     X_s, loc_slices, domain_priors,
                                     y, b, gp, gp_prior,
                                     prior_jump_probability,
                                     log_alt = bmim_joint_log_alt_kernel) {
  if (xi[p] == 0L || b == 0L) {
    # The data are uninformative about A_p when xi[p]=0 or b=0.
    # Draw A_p exactly from its prior. IMPORTANT: if b=0 but xi[p]=1,
    # the ACTIVE distance sum must still be refreshed, because the next
    # Bernoulli update evaluates the GP alternative using this same A_p.
    new_A <- bmim_domain_draw(domain_priors[[p]])
    new_Dp <- bmim_domain_distance(X_s, loc_slices[[p]], new_A)
    if (xi[p] == 1L) {
      D_total <- bmim_joint_distance_replace(
        D_total, D_by_loc[[p]], new_Dp
      )
    }
    A[p, ] <- new_A
    D_by_loc[[p]] <- new_Dp
    return(list(A = A, D_by_loc = D_by_loc,
                D_total = D_total, accepted = TRUE))
  }
  proposed <- bmim_domain_propose(A[p, ], domain_priors[[p]],
                                   prior_jump_probability)
  if (identical(proposed$kind, "identity") ||
      all(proposed$A == A[p, ])) {
    return(list(A = A, D_by_loc = D_by_loc,
                D_total = D_total, accepted = FALSE))
  }
  D_new_p <- bmim_domain_distance(X_s, loc_slices[[p]], proposed$A)
  D_proposed <- bmim_joint_distance_replace(
    D_total, D_by_loc[[p]], D_new_p
  )
  lr <- log_alt(y, D_proposed, gp, gp_prior) -
    log_alt(y, D_total, gp, gp_prior) +
    proposed$log_prior_and_proposal_ratio
  accepted <- is.finite(lr) && log(runif(1L)) < min(0, lr)
  if (accepted) {
    A[p, ] <- proposed$A
    D_by_loc[[p]] <- D_new_p
    D_total <- D_proposed
  }
  list(A = A, D_by_loc = D_by_loc,
       D_total = D_total, accepted = accepted)
}

# Joint node/interval proposal. Propose xi_p* = 1-xi_p and simultaneously
# draw A_p* from its prior. Because q(A_p*) = p(A_p*) in BOTH directions,
# the domain prior cancels the reverse/forward proposal correction EXACTLY.
# The remaining MH ratio is the GP likelihood plus Ising/VB prior ratio.
bmim_joint_flip_node <- function(p, xi, A, D_by_loc, D_total,
                                 X_s, loc_slices, domain_priors,
                                 y, b, gp, gp_prior, W_mat, p_prop,
                                 log_alt = bmim_joint_log_alt_kernel) {
  xi_new <- xi
  xi_new[p] <- 1L - xi[p]
  if (sum(xi_new) == 0L) {
    return(list(xi = xi, A = A, D_by_loc = D_by_loc,
                D_total = D_total, accepted = FALSE, attempted = FALSE))
  }
  A_new_p <- bmim_domain_draw(domain_priors[[p]])
  D_new_p <- bmim_domain_distance(X_s, loc_slices[[p]], A_new_p)
  D_prop <- D_total
  if (xi[p] == 1L) {
    D_prop <- bmim_joint_distance_replace(D_prop, D_by_loc[[p]],
                                           matrix(0, nrow(D_prop), ncol(D_prop)))
  }
  if (xi_new[p] == 1L) {
    D_prop <- bmim_joint_distance_replace(D_prop,
                                           matrix(0, nrow(D_prop), ncol(D_prop)),
                                           D_new_p)
  }
  lr <- bmim_joint_log_xi_prior(xi_new, W_mat, p_prop) -
    bmim_joint_log_xi_prior(xi, W_mat, p_prop)
  if (b == 1L) {
    lr <- lr + log_alt(y, D_prop, gp, gp_prior) -
      log_alt(y, D_total, gp, gp_prior)
  }
  accepted <- is.finite(lr) && log(runif(1L)) < min(0, lr)
  if (accepted) {
    xi <- xi_new
    A[p, ] <- A_new_p
    D_by_loc[[p]] <- D_new_p
    D_total <- D_prop
  }
  list(xi = xi, A = A, D_by_loc = D_by_loc,
       D_total = D_total, accepted = accepted, attempted = TRUE)
}

