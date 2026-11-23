library(dplyr)
library(tidyr)

node <- read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_node_results_30.csv")
domain <- read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_domain_results_30.csv")

# =====================================================
# Helper: compute metrics
# =====================================================

metrics <- function(TP,FP,FN,TN){
  
  TPR <- TP/(TP+FN)
  FPR <- FP/(FP+TN)
  PPV <- TP/(TP+FP)
  F1  <- 2*TP/(2*TP+FP+FN)
  
  c(TPR=TPR,FPR=FPR,PPV=PPV,F1=F1)
}

# =====================================================
# Aggregate node confusion matrices
# =====================================================

node_summary <- node %>%
  group_by(method) %>%
  summarise(
    TP=sum(TP,na.rm=TRUE),
    FP=sum(FP,na.rm=TRUE),
    FN=sum(FN,na.rm=TRUE),
    .groups="drop"
  )

num_reps <- length(unique(node$replicate))
num_locs <- 30

total <- num_locs*num_reps

node_summary <- node_summary %>%
  mutate(TN = total - TP - FP - FN)

node_summary <- node_summary %>%
  rowwise() %>%
  mutate(vals=list(metrics(TP,FP,FN,TN))) %>%
  unnest_wider(vals)

# =====================================================
# Identify best values
# =====================================================

best_TPR <- max(node_summary$TPR,na.rm=TRUE)
best_FPR <- min(node_summary$FPR,na.rm=TRUE)
best_PPV <- max(node_summary$PPV,na.rm=TRUE)
best_F1  <- max(node_summary$F1,na.rm=TRUE)

fmt <- function(x) sprintf("%.3f",x)

latex_val <- function(x,best,maximize=TRUE){
  
  if(is.na(x)) return("-")
  
  if(maximize){
    if(abs(x-best)<1e-10) return(paste0("\\bf ",fmt(x)))
  } else {
    if(abs(x-best)<1e-10) return(paste0("\\bf ",fmt(x)))
  }
  
  fmt(x)
}

# =====================================================
# Choose the two BMIM node variants
# (modify names if needed)
# =====================================================

bmim1 <- node_summary[node_summary$method=="A3_Domain_ISGM",]
bmim2 <- node_summary[node_summary$method=="A4_ISGM_Full",]

# =====================================================
# Print LaTeX rows
# =====================================================

cat("\nPASTE THIS INTO OVERLEAF:\n\n")

cat(
  "$\\mathbf{TPR}$ &",
  latex_val(bmim1$TPR,best_TPR,TRUE),"&",
  latex_val(bmim2$TPR,best_TPR,TRUE),"\\\\\n"
)

cat(
  "$\\mathbf{FPR}$ &",
  latex_val(bmim1$FPR,best_FPR,FALSE),"&",
  latex_val(bmim2$FPR,best_FPR,FALSE),"\\\\\n"
)

cat(
  "$\\mathbf{PPV}$ &",
  latex_val(bmim1$PPV,best_PPV,TRUE),"&",
  latex_val(bmim2$PPV,best_PPV,TRUE),"\\\\\n"
)

cat(
  "$\\mathbf{F_{1}}$ &",
  latex_val(bmim1$F1,best_F1,TRUE),"&",
  latex_val(bmim2$F1,best_F1,TRUE),"\\\\\n"
)