# =========================================================
# correlation_matrix_constructor.R
# Construction of W matrices for ISGM
# =========================================================

# ---------------------------------------------------------
# Pairwise functional correlation between two locations
# ---------------------------------------------------------
pairwise_functional_cor <- function(mat_l, mat_k, method = "pearson",
                                    agg = c("mean", "median")){
  agg <- match.arg(agg)
  
  stopifnot(
    nrow(mat_l) == nrow(mat_k),
    ncol(mat_l) == ncol(mat_k)
  )
  
  N <- nrow(mat_l)
  cors <- numeric(N)
  
  for(i in seq_len(N)){
    cors[i] <- suppressWarnings(
      cor(mat_l[i, ], mat_k[i, ], method = method, use = "complete.obs")
    )
  }
  
  if(agg == "mean"){
    mean(cors, na.rm = TRUE)
  } else {
    median(cors, na.rm = TRUE)
  }
}




#' Build W matrices for ISGM
#'
#' @param mode "simulation" or "data"
#' @param full_dat2 list of (N x T) data.frames, one per location
#' @param chosen_locs integer vector of true signal locations (simulation only)
#' @param t_length number of time points
#' @param scale_full scaling factor for W_func_full
#'
#' @return list(W_func, W_func_full)
#' @export
build_W_matrix <- function(
    mode = c("simulation", "data"),
    full_dat2,
    chosen_locs = NULL,
    t_length,
    scale_full = 0.4,
    cor_method = "pearson",
    agg = c("mean", "median")
){
  
  agg <- match.arg(agg)
  mode <- match.arg(mode)
  num_locs <- length(full_dat2)
  
  W_func <- matrix(0, nrow = num_locs, ncol = num_locs)
  W_func_full <- matrix(0, nrow = num_locs, ncol = num_locs)
  
  if(mode == "simulation"){
    if(is.null(chosen_locs)){
      stop("chosen_locs must be provided in simulation mode")
    }
    
    for(l in seq_len(num_locs)){
      for(k in seq_len(num_locs)){
        if(l >= k) next
        
        if((l %in% chosen_locs) && (k %in% chosen_locs)){
          # 🔥 artificially strong correlation
          W_func[k,l] <- rnorm(1, 0.75, 0.05)
          W_func_full[k,l] <- scale_full * W_func[k,l]
        } else {
          # weak / noisy correlation
          tmp <- rt(1, t_length - 2, 0)
          tmp <- tmp * ((1 - tmp^2) / (t_length - 2))
          W_func[k,l] <- tmp
          W_func_full[k,l] <- tmp
        }
        
        W_func[l,k] <- W_func[k,l]
        W_func_full[l,k] <- W_func_full[k,l]
      }
    }
  }
  
  if(mode == "data"){
    for(l in seq_len(num_locs)){
      for(k in seq_len(num_locs)){
        if(l >= k) next
        
        mat_l <- as.matrix(full_dat2[[l]])
        mat_k <- as.matrix(full_dat2[[k]])
        
        rho_lk <- pairwise_functional_cor(
          mat_l, mat_k,
          method = cor_method,
          agg = agg
        )
        
        W_func[k,l] <- rho_lk
        W_func_full[k,l] <- scale_full * rho_lk
        
        W_func[l,k] <- W_func[k,l]
        W_func_full[l,k] <- W_func_full[k,l]
      }
    }
  }
  
  list(
    W_func = W_func,
    W_func_full = W_func_full
  )
}