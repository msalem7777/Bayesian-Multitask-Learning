# =========================================================
# location_selection.R
# Spatial/location selection using ISGM
# =========================================================

# Requires:
#   ISGM() implementation available in your environment
#   (typically from your Ising selection code)


#' Run ISGM-based location selection
#'
#' @param final_data list of data.frames (after domain truncation), per location
#' @param sim_scores numeric response vector (length N)
#' @param W_func adjacency / similarity matrix between locations (L x L)
#' @param prior_adj scalar prior strength multiplier
#' @param grps_vec grouping vectors from domain selection
#' @param num_locs number of locations
#' @param N number of MCMC / VB iterations
#' @param VBprop proposal type flag
#'
#' @return xi_chain matrix (iterations x num_locs)
#'
#' @export
run_location_selection <- function(final_data,
                                   sim_scores,
                                   W_func,
                                   prior_adj = 1,
                                   grps_vec,
                                   num_locs,
                                   N = 300,
                                   VBprop = 1){

  # ---- build integrated design matrix ----
  # rows: individuals, columns: stacked (location × time)

  final_int_df <- NULL
  for(l in seq_len(num_locs)){
    final_int_df <- cbind(final_int_df, final_data[[l]])
  }
  final_int_df <- as.data.frame(final_int_df)

  # ---- initial xi ----
  curr_xi <- rep(1, num_locs)

  # ---- scale adjacency ----
  W_mat <- prior_adj * W_func

  # ---- ISGM call ----
  xi_chain <- ISGM(
    path_features = final_int_df,
    W_mat = W_mat,
    y = sim_scores,
    curr_xi = curr_xi,
    grps_vec = grps_vec,
    connec = 1.0,
    N = N,
    VBprop = VBprop
  )

  xi_chain
}


# =========================================================
# Baseline: no domain selection (full features)
# =========================================================

#' Run ISGM without domain truncation
#' @export
run_location_selection_full <- function(full_dat2,
                                        sim_scores,
                                        W_func,
                                        prior_adj = 1,
                                        num_locs,
                                        N = 300,
                                        VBprop = 1){

  final_int_df <- NULL
  for(l in seq_len(num_locs)){
    final_int_df <- cbind(final_int_df, full_dat2[[l]])
  }
  final_int_df <- as.data.frame(final_int_df)

  curr_xi <- rep(1, num_locs)
  W_mat <- prior_adj * W_func

  grps_full <- build_full_grps(full_dat2)

  xi_chain <- ISGM(
    path_features = final_int_df,
    W_mat = W_mat,
    y = sim_scores,
    curr_xi = curr_xi,
    grps_vec = grps_full,
    connec = 1.0,
    N = N,
    VBprop = VBprop
  )

  xi_chain
}
