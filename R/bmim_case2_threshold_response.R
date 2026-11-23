# Case 2 mathematical response rule from the manuscript (eq:mod29).
# This is a PURE response generator once K_gen, alpha, and empirical moments
# have been supplied. Construction of K_gen and alpha is NOT supplied by the
# existing historical R simulation source. Do NOT substitute an arbitrary
# kernel or coefficient distribution and present it as manuscript Case 2.
#
# Positive side: eta_i = [K_gen %*% alpha]_i >= 0 => group 1.
# Negative side: eta_i < 0 => group 2.
# Null: caller must generate a response independent of X instead; this
# function deliberately implements only the stated GP-alternative rule.
bmim_case2_threshold_response <- function(K_gen, alpha,
                                          mu_positive, sd_positive,
                                          mu_negative, sd_negative) {
  if (!is.matrix(K_gen) || !is.numeric(K_gen) ||
      nrow(K_gen) < 2L || ncol(K_gen) < 1L ||
      any(!is.finite(K_gen)) ||
      !is.numeric(alpha) || length(alpha) != ncol(K_gen) ||
      any(!is.finite(alpha))) {
    stop("K_gen must be a finite numeric matrix and alpha a compatible numeric vector.")
  }
  moments <- c(mu_positive, sd_positive, mu_negative, sd_negative)
  if (length(moments) != 4L || any(!is.finite(moments)) ||
      sd_positive <= 0 || sd_negative <= 0) {
    stop("Two finite group means and strictly positive group SDs are required.")
  }
  score <- drop(K_gen %*% alpha)
  group_positive <- score >= 0
  mu <- ifelse(group_positive, mu_positive, mu_negative)
  sd <- ifelse(group_positive, sd_positive, sd_negative)
  list(y = rnorm(nrow(K_gen), mean = mu, sd = sd),
       score = score, group_positive = group_positive,
       mean = mu, sd = sd)
}
