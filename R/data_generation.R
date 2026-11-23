# =========================================================
# data_generation.R
# Functional data simulation utilities
# =========================================================
library(mvtnorm)
#' Create base functional trajectories per location
#' @export
fn_dat_creator <- function(num_locs, num_ind, num_T, dist_mat, gt_tau){
  listr <- vector("list", num_ind)
  for(i in seq_len(num_ind)){
    Wiener <- matrix(NA, nrow=num_locs, ncol=num_T)
    for(j in seq_len(num_T)){
      Wiener[,j] <- mvtnorm::rmvnorm(
        1,
        mean = rep(0, num_locs),
        sigma = gt_tau * exp(-0.5*(dist_mat - mean(dist_mat)) / sd(dist_mat))
      )
    }
    # Wiener_dt = apply(Wiener, 1, detrend)
    listr[[i]] <- Wiener
  }
  return(listr)
}

#' Append irrelevant functional segment
#' @export
fn_dat_appendor <- function(num_locs, num_ind, num_T, og_dfs, dist_mat, gt_sigma_u){
  listr <- vector("list", num_ind)
  for(i in seq_len(num_ind)){
    Wiener <- matrix(NA, nrow=num_locs, ncol=num_T)
    for(j in seq_len(num_T)){
      if(j == 1){
        Wiener[,j] <- mvtnorm::rmvnorm(
          1, rep(0, num_locs),
          gt_sigma_u * exp(-0.5*(dist_mat - mean(dist_mat)) / sd(dist_mat))
        ) + og_dfs[[i]][,num_T]
      } else{
        Wiener[,j] <- mvtnorm::rmvnorm(
          1, rep(0, num_locs),
          gt_sigma_u * exp(-0.5*(dist_mat - mean(dist_mat)) / sd(dist_mat))
        )
      }
    }
    listr[[i]] <- Wiener
  }
  return(listr)
}

#' Prepend irrelevant functional segment
#' @export
fn_dat_preppendor <- function(num_locs, num_ind, num_T, og_dfs, dist_mat, gt_sigma_u){
  listr <- vector("list", num_ind)
  for(i in seq_len(num_ind)){
    Wiener <- matrix(NA, nrow=num_locs, ncol=num_T)
    for(j in num_T:1){
      if(j == num_T){
        Wiener[,j] <- og_dfs[[i]][,1] - mvtnorm::rmvnorm(
          1, rep(0, num_locs),
          gt_sigma_u * exp(-0.5*(dist_mat - mean(dist_mat)) / sd(dist_mat))
        )
      } else{
        Wiener[,j] <- - mvtnorm::rmvnorm(
          1, rep(0, num_locs),
          gt_sigma_u * exp(-0.5*(dist_mat - mean(dist_mat)) / sd(dist_mat))
        )
      }
    }
    listr[[i]] <- Wiener
  }
  return(listr)
}

#' Full dataset simulator (all stages)
#' @export
simulate_dataset <- function(config, real_process, chosen_locs){

  with(config, {

    chosen_one   <- 1
    unchosen_one <- 2

    true_locs <- rep(0, num_locs)
    if(real_process == 1) true_locs[chosen_locs] <- 1

    # ---- distance matrix ----
    dist_mat <- matrix(0, num_locs, num_locs)
    for(i in seq_len(num_locs)){
      for(j in i:num_locs){
        if(i == j) next
        if(i %in% chosen_locs && j %in% chosen_locs){
          dist_mat[i,j] <- dist_mat[j,i] <- rgamma(1, 4.5, 3.5)
        } else{
          dist_mat[i,j] <- dist_mat[j,i] <- rnorm(1, 4.5, 1.5)
        }
      }
    }

    # ---- base functional data ----
    fn_dat_dfs <- fn_dat_creator(num_locs, 2, num_T, dist_mat, gt_tau)
    for(i in 3:nx) fn_dat_dfs[[i]] <- matrix(0, num_locs, num_T)

    true_status <- numeric(nx)
    true_status[1:2] <- c(1,0)

    for(i in 3:nx){
      true_status[i] <- rbinom(1,1,0.5)
      for(j in seq_len(num_locs)){
        if(real_process == 1 && true_status[i] == 1 && j %in% chosen_locs){
          fn_dat_dfs[[i]][j,] <- fn_dat_dfs[[chosen_one]][j,] + rnorm(num_T, 0, gt_sigma_pop)
        } else{
          fn_dat_dfs[[i]][j,] <- fn_dat_dfs[[unchosen_one]][j,] + rnorm(num_T, 0, gt_sigma_pop)
        }
      }
    }

    # ---- append + prepend ----
    append1 <- fn_dat_appendor(num_locs, 2, num_T, fn_dat_dfs, dist_mat, gt_sigma_u)
    prep1   <- fn_dat_preppendor(num_locs, 2, num_T, fn_dat_dfs, dist_mat, gt_sigma_u)

    for(i in 3:nx){
      append1[[i]] <- matrix(0, num_locs, num_T)
      prep1[[i]]   <- matrix(0, num_locs, num_T)
    }

    # ---- combine stages ----
    full_dat <- Map(cbind, fn_dat_dfs, append1)
    full_dat <- Map(cbind, full_dat, prep1)

    # ---- per-location wide format ----
    full_dat2 <- vector("list", num_locs)
    for(l in seq_len(num_locs)){
      mat <- matrix(0, nx, ncol(full_dat[[1]]))
      for(i in seq_len(nx)) mat[i,] <- full_dat[[i]][l,]
      full_dat2[[l]] <- as.data.frame(mat)
    }

    list(
      full_dat2   = full_dat2,
      true_status = true_status,
      true_locs   = true_locs,
      dist_mat    = dist_mat
    )
  })
}
