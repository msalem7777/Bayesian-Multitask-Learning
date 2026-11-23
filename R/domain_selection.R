# =========================================================
# domain_selection.R
# Time-domain selection using FocusGP
# =========================================================

# Requires:
#   FocusGP(), FCV(), EI.search(), etc. (from your GP optimization code)
#   These should live in: domain_selection_core.R or gp_utils.R


#' Run time-domain selection for each location
#'
#' @param full_dat2 list of data.frames, one per location (N x T)
#' @param sim_scores numeric response (length N)
#' @param t_length number of time points per location
#' @param num_locs number of locations
#'
#' @return list with:
#'   final_data: truncated + reordered predictors per location
#'   grps_vec: list of consecutive index groups per location
#'   domains: binary vectors of selected time indices per location
#'
#' @export
#' Run time-domain selection for each location (original-faithful)
#'
#' Mimics the original pipeline:
#'  - Runs FocusGP per location to pick [a_sel, b_sel]
#'  - Centers each subject trajectory (scale=FALSE) over ALL t_length columns
#'  - Keeps ONLY the selected window (not reordering leftovers)
#'  - Appends an intercept column = per-subject mean over ALL t_length
#'  - Stores grps_vec as consecutive runs of selected indices
#'  - Stores domains as 0/1 vectors length t_length
#'
#' @export
run_domain_selection <- function(full_dat2, sim_scores, t_length, num_locs,
                                 nfold = 5, lambda = 0.05, ninit = 12, end = 30,
                                 max_retries = 5){
  
  stopifnot(length(full_dat2) == num_locs)
  N <- length(sim_scores)
  
  # outputs
  grps_vec <- vector("list", num_locs)
  dat_listr_domains <- vector("list", num_locs)
  final_data <- vector("list", num_locs)
  
  for(l in seq_len(num_locs)){
    
    X <- as.matrix(full_dat2[[l]][, 1:t_length, drop = FALSE])
    Z <- sim_scores
    
    # ---- FocusGP with retry (no deleteGPsep(gpi) here; FocusGP manages its own gpi)
    out <- NULL
    for(r in seq_len(max_retries)){
      out <- tryCatch(
        FocusGP(X, Z, nfold = nfold, lambda = lambda, ninit = ninit, end = end),
        error = function(e){
          message("FocusGP failed at location ", l, " (retry ", r, "): ", e$message)
          NULL
        }
      )
      if(!is.null(out)) break
    }
    if(is.null(out)){
      stop("FocusGP failed after ", max_retries, " retries at location ", l)
    }
    
    # ---- selected window
    idx <- which.min(out$maxei)
    
    a_sel <- out$X[idx, 1]
    b_sel <- a_sel + out$X[idx, 2]
    
    # original code didn't clamp explicitly; do it to avoid OOB
    a_sel <- max(1, floor(a_sel))
    b_sel <- min(t_length, floor(b_sel))
    
    GP_select <- seq(a_sel, b_sel)
    
    # ---- group consecutive indices
    consec_vecs <- split(GP_select, cumsum(c(1, diff(GP_select) != 1)))
    grps_vec[[l]] <- consec_vecs
    
    # ---- domain indicator
    dom <- rep(0, t_length)
    dom[GP_select] <- 1
    dat_listr_domains[[l]] <- dom
    
    # ---- center each subject across ALL t_length timepoints (scale=FALSE)
    appnder <- as.data.frame(t(apply(full_dat2[[l]][, 1:t_length, drop = FALSE], 1, scale, scale = FALSE)))
    colnames(appnder) <- paste0("t", seq_len(t_length))
    
    # ---- keep ONLY selected window
    appnder_sel <- appnder[, GP_select, drop = FALSE]
    
    # ---- intercepts (per-subject mean over ALL t_length, using original logic)
    intrcpts <- rowMeans(full_dat2[[l]][, 1:t_length, drop = FALSE], na.rm = TRUE)
    
    # ---- final truncated design for this location = selected window + intercept
    final_data[[l]] <- cbind(as.data.frame(appnder_sel), intrcpts = intrcpts)
  }
  
  # integrated frames (like original)
  final_int_df <- as.data.frame(final_data)
  
  # also build "full" branch like your current function (for ISGM full)
  final_data_full <- lapply(full_dat2, function(df) as.data.frame(df)[, 1:t_length, drop = FALSE])
  grps_vec_full <- lapply(final_data_full, function(df) list(seq_len(ncol(df))))
  final_int_df_full <- as.data.frame(final_data_full)
  
  list(
    final_data = final_data,
    final_int_df = final_int_df,
    grps_vec = grps_vec,
    domains = dat_listr_domains,
    final_data_full = final_data_full,
    final_int_df_full = final_int_df_full,
    grps_vec_full = grps_vec_full
  )
}



# =========================================================
# Optional helper: build full-domain grouping (no truncation)
# =========================================================

#' Build grouping vectors when no domain selection is used
#' @export
build_full_grps <- function(full_dat2){
  lapply(full_dat2, function(df) list(seq_len(ncol(df))))
}
