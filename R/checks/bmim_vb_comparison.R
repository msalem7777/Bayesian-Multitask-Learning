# Optional paired engineering pilot, NOT the paper's production experiment.
# Run from R/: Rscript --vanilla -e "source('checks/bmim_vb_comparison.R')"
# Set BMIM_VB_PILOT_REPS and BMIM_VB_PILOT_ITER for larger comparisons.
# Uses identical data, network, GP priors and chain seeds across measures.
files <- Filter(Negate(is.null),lapply(sys.frames(),function(f) f$ofile))
root <- dirname(dirname(normalizePath(tail(files,1L)[[1L]])))
source(file.path(root,"bmim_gp_updates.R"))
source(file.path(root,"bmim_global_gp_test.R"))
source(file.path(root,"bmim_joint_domains.R"))
for (pkg in c("ClusterR","infotheo"))
  if (!requireNamespace(pkg,quietly=TRUE)) stop("Pilot requires ",pkg)
read_int <- function(key,default,lower) {
  x <- suppressWarnings(as.integer(Sys.getenv(key,unset=as.character(default))))
  if (is.na(x) || x<lower) stop("Invalid ",key)
  x
}
reps <- read_int("BMIM_VB_PILOT_REPS",2L,1L)
iters <- read_int("BMIM_VB_PILOT_ITER",300L,20L)
# Small illustrative problem: 12 nodes, 6 measurements, 2 active nodes.
# With Z, node 3 is a confounded proxy for Z; the actual signal is nonlinear.
# Each scenario uses a new dataset; methods within a scenario share it.
rows <- list()
for (rep in seq_len(reps)) for (has_z in c(FALSE,TRUE)) for (signal in c(FALSE,TRUE)) {
  set.seed(70000L+100L*rep+10L*as.integer(has_z)+as.integer(signal))
  n <- 100L; P <- 12L; T <- 6L
  z <- rnorm(n); latent <- matrix(rnorm(n*P),n,P)
  blocks <- lapply(seq_len(P),function(p) {
    x <- matrix(rnorm(n*T),n,T)
    x[,2:4] <- x[,2:4]+latent[,p]
    if (has_z && p==3L) x <- x+2*z
    x
  })
  X <- do.call(cbind,blocks)
  loc <- lapply(seq_len(P),function(p) seq.int((p-1L)*T+1L,p*T))
  groups <- rep(list(seq_len(T)),P)
  u <- rowMeans(blocks[[1L]][,2:4]); v <- rowMeans(blocks[[2L]][,2:4])
  y <- rnorm(n,sd=.7)+(if (has_z) 2*z else 0)
  if (signal) y <- y+as.numeric(scale(u^2))+as.numeric(scale(sin(v)))
  Z <- if (has_z) matrix(z,ncol=1) else NULL
  W <- matrix(0,P,P) # neutral network isolates the effect of association measure
  methods <- if (has_z) c("sq.corr","partial.sq.corr","mi","residual.mi") else c("sq.corr","mi")
  true_nodes <- if (signal) as.integer(seq_len(P)<=2L) else integer(P)
  true_time <- matrix(FALSE,P,T); if (signal) true_time[1:2,2:4] <- TRUE
  for (method in methods) {
    start <- proc.time()[["elapsed"]]
    vb <- bmim_vb_prior(X,y,loc,W,Z,method,seed=71000L+rep)
    prep <- proc.time()[["elapsed"]]-start
    start <- proc.time()[["elapsed"]]
    fit <- bmim_joint_domains_from_features(X,W,y,groups,p_prop=vb$weights,
      max_domain_length=3L,Z=Z,n_iter=iters,burnin=iters%/%4L,
      n_chains=4L,pi_b=.5,seed=72000L+rep)
    runtime <- proc.time()[["elapsed"]]-start
    # Joint probabilities account for global absence; conditional probabilities
    # would include only iterations in which the functional GP is present.
    pnode <- fit$posterior_gp_and_node_included
    selected <- pnode>.5; selected_time <- fit$posterior_node_time_joint_gp>.5
    tp <- sum(selected & true_nodes==1L); fp <- sum(selected & true_nodes==0L)
    union <- sum(selected_time | true_time)
    switches <- sum(vapply(fit$chains,function(ch) {
      b <- ch$parameters[seq.int(fit$burnin+1L,iters),"b"]
      sum(diff(b)!=0)
    },numeric(1L)))
    rows[[length(rows)+1L]] <- data.frame(rep=rep,scalar=has_z,signal=signal,method=method,
      gp_probability=fit$posterior_gp_included,global_brier=(fit$posterior_gp_included-signal)^2,
      node_brier=mean((pnode-true_nodes)^2),tpr=if(signal) tp/2 else NA_real_,
      fpr=fp/sum(true_nodes==0L),ppv=if(sum(selected)) tp/sum(selected) else NA_real_,
      f1=if(signal) 2*tp/(sum(selected)+2) else NA_real_,
      iou=if(union) sum(selected_time & true_time)/union else NA_real_,
      active_weight=if(signal) mean(vb$weights[1:2]) else NA_real_,
      inactive_weight=mean(vb$weights[true_nodes==0L]),global_switches=switches,
      vb_seconds=prep,mcmc_seconds=runtime)
    print(tail(rows,1)[[1L]],row.names=FALSE)
  }
}
results <- do.call(rbind,rows)
# No output files unless explicitly requested. Avoid overwriting old results.
path <- Sys.getenv("BMIM_VB_PILOT_OUTPUT",unset="")
if(nzchar(path)) {
  if(file.exists(path)) stop("Pilot output exists; select a new path.")
  write.csv(results,path,row.names=FALSE)
}
cat("BMIM VB PAIRED PILOT COMPLETE; exploratory only, not a calibration study.\n")
