# =========================================================
# evaluation_metrics.R
# =========================================================

#' @export
selection_metrics <- function(selected, truth) {
  selected <- as.logical(selected)
  truth    <- as.logical(truth)
  
  TP <- sum(selected & truth)
  FP <- sum(selected & !truth)
  FN <- sum(!selected & truth)
  
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), NA)
  recall    <- ifelse(TP + FN > 0, TP / (TP + FN), NA)
  F1        <- ifelse(is.na(precision + recall), NA,
                      2 * precision * recall / (precision + recall))
  
  c(
    n_selected = sum(selected),
    TP = TP,
    FP = FP,
    FN = FN,
    precision = precision,
    recall = recall,
    F1 = F1
  )
}

#' @export
log_method <- function(method, selected_locs, tau, truth, meta, rep_id) {
  
  level <- "Location"
  
  if (all(is.na(selected_locs)) || all(is.na(truth))) {
    TP <- FP <- FN <- precision <- recall <- F1 <- NA_real_
    n_selected <- if (all(is.na(selected_locs))) NA_integer_ else sum(as.logical(selected_locs))
  } else {
    selected_locs <- as.logical(selected_locs)
    truth <- as.logical(truth)
    
    TP <- sum(selected_locs & truth)
    FP <- sum(selected_locs & !truth)
    FN <- sum(!selected_locs & truth)
    
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    F1 <- if (is.finite(precision + recall) && (precision + recall) > 0)
      2 * precision * recall / (precision + recall)
    else NA_real_
    
    n_selected <- sum(selected_locs)
  }
  
  data.frame(
    replicate    = rep_id,
    level        = level,
    method       = method,
    tau          = tau,
    n_selected   = n_selected,
    TP           = TP,
    FP           = FP,
    FN           = FN,
    precision    = precision,
    recall       = recall,
    F1           = F1,
    real_process = meta$real_process,
    stringsAsFactors = FALSE
  )
}

#' @export
log_domain_method <- function(
    method,
    selected_domains,   # matrix [L × T]
    truth_domains,      # matrix [L × T]
    meta,
    rep_id
) {
  stopifnot(
    is.matrix(selected_domains),
    is.matrix(truth_domains),
    all(dim(selected_domains) == dim(truth_domains))
  )
  
  sel <- as.logical(selected_domains)
  tru <- as.logical(truth_domains)
  
  TP <- sum(sel & tru)
  FP <- sum(sel & !tru)
  FN <- sum(!sel & tru)
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  F1 <- if (!is.na(precision + recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall)
  else NA_real_
  
  data.frame(
    replicate    = rep_id,
    level        = "Domain",
    method       = method,
    tau          = NA_real_,
    n_selected   = sum(sel),
    TP           = TP,
    FP           = FP,
    FN           = FN,
    precision    = precision,
    recall       = recall,
    F1           = F1,
    real_process = meta$real_process,
    stringsAsFactors = FALSE
  )
}

#' @export
iou_binary <- function(pred, truth) {
  pred  <- as.logical(pred)
  truth <- as.logical(truth)
  
  inter <- sum(pred & truth)
  union <- sum(pred | truth)
  
  if (union == 0) return(NA_real_)
  inter / union
}

#' @export
log_domain_iou <- function(method, pred_domains, truth_time_mask, meta, rep_id) {
  
  # pred_domains: list of length num_locs
  # each element is integer vector of selected time indices (or NULL)
  
  ious <- numeric(length(pred_domains))
  
  for (l in seq_along(pred_domains)) {
    pred_mask <- rep(FALSE, length(truth_time_mask))
    
    if (!is.null(pred_domains[[l]]) && length(pred_domains[[l]]) > 0) {
      pred_mask[pred_domains[[l]]] <- TRUE
    }
    
    ious[l] <- iou_binary(pred_mask, truth_time_mask)
  }
  
  data.frame(
    replicate    = rep_id,
    level        = "Domain",
    method       = method,
    mean_IoU     = mean(ious, na.rm = TRUE),
    median_IoU   = median(ious, na.rm = TRUE),
    n_loc_active = sum(!is.na(ious)),
    real_process = meta$real_process,
    stringsAsFactors = FALSE
  )
}


# collapse time-domain mask to location-level decision
collapse_domain <- function(selected_mask, true_mask) {
  sel_loc  <- vapply(seq_along(selected_mask),
                     function(j) any(selected_mask[[j]]), logical(1))
  true_loc <- vapply(seq_along(true_mask),
                     function(j) any(true_mask[[j]]), logical(1))
  
  TP <- sum(sel_loc & true_loc)
  FP <- sum(sel_loc & !true_loc)
  FN <- sum(!sel_loc & true_loc)
  TN <- sum(!sel_loc & !true_loc)
  
  c(TP=TP, FP=FP, FN=FN, TN=TN)
}

# node confusion
node_confusion <- function(selected, truth) {
  TP <- sum(selected & truth)
  FP <- sum(selected & !truth)
  FN <- sum(!selected & truth)
  TN <- sum(!selected & !truth)
  c(TP=TP, FP=FP, FN=FN, TN=TN)
}


# -------------------------------
# IoU metrics for time-domain masks
# -------------------------------
#' @export
domain_metrics_iou <- function(selected, truth) {
  selected <- as.logical(selected)
  truth    <- as.logical(truth)
  
  stopifnot(length(selected) == length(truth))
  
  TP <- sum(selected & truth)
  FP <- sum(selected & !truth)
  FN <- sum(!selected & truth)
  denom <- TP + FP + FN
  
  iou <- if (denom > 0) TP / denom else NA_real_
  
  c(
    n_selected = sum(selected),
    TP = TP,
    FP = FP,
    FN = FN,
    IoU = iou
  )
}

# -------------------------------
# Log domain results PER LOCATION
# selected_domains_by_loc: list length num_locs, each logical length t_length
# truth_domains_by_loc:    list length num_locs, each logical length t_length
# -------------------------------
#' @export
log_domain_method_per_loc <- function(method, selected_domains_by_loc, truth_domains_by_loc, meta, rep_id) {
  
  L <- length(truth_domains_by_loc)
  stopifnot(length(selected_domains_by_loc) == L)
  
  rows <- vector("list", L)
  
  for (l in seq_len(L)) {
    sel <- selected_domains_by_loc[[l]]
    tru <- truth_domains_by_loc[[l]]
    
    # allow NA for methods with no domain selection, or no ground truth
    if (all(is.na(sel)) || all(is.na(tru))) {
      rows[[l]] <- data.frame(
        replicate = rep_id,
        level = "Domain",
        method = method,
        location = l,
        n_selected = if (all(is.na(sel))) NA_integer_ else sum(as.logical(sel)),
        TP = NA_integer_,
        FP = NA_integer_,
        FN = NA_integer_,
        IoU = NA_real_,
        real_process = meta$real_process,
        stringsAsFactors = FALSE
      )
    } else {
      m <- domain_metrics_iou(sel, tru)
      rows[[l]] <- data.frame(
        replicate = rep_id,
        level = "Domain",
        method = method,
        location = l,
        n_selected = unname(m["n_selected"]),
        TP = unname(m["TP"]),
        FP = unname(m["FP"]),
        FN = unname(m["FN"]),
        IoU = unname(m["IoU"]),
        real_process = meta$real_process,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, rows)
}

# -------------------------------
# Build truth domain masks per location
# Assumption (matches your note): true time window is indices 61:80,
# and only chosen_locs are "truly active" locations when real_process == 1.
# -------------------------------
#' @export
make_truth_domains_by_loc <- function(num_locs, t_length, chosen_locs, real_process,
                                      true_idx = 61:80) {
  out <- vector("list", num_locs)
  for (l in seq_len(num_locs)) {
    mask <- rep(FALSE, t_length)
    if (real_process == 1 && (l %in% chosen_locs)) {
      mask[true_idx] <- TRUE
    }
    out[[l]] <- mask
  }
  out
}

