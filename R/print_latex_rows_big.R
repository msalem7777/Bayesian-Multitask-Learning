library(dplyr)
library(tidyr)

# =========================================================
# LOAD RESULTS
# =========================================================

node <- read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_node_results_30.csv")
domain <- read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_domain_results_30.csv")

# =========================================================
# HELPERS
# =========================================================

fmt <- function(x) sprintf("%.3f", x)

bold_max <- function(x, best){
  if(is.na(x)) return("-")
  if(abs(x-best) < 1e-10) return(paste0("\\bf ", fmt(x)))
  fmt(x)
}

bold_min <- function(x, best){
  if(is.na(x)) return("-")
  if(abs(x-best) < 1e-10) return(paste0("\\bf ", fmt(x)))
  fmt(x)
}

safe_val <- function(df, method, col){
  v <- df[df$method == method, col]
  if(length(v)==0) return(NA)
  v
}

metrics <- function(TP,FP,FN,TN){
  TPR <- TP/(TP+FN)
  FPR <- FP/(FP+TN)
  PPV <- TP/(TP+FP)
  F1  <- 2*TP/(2*TP+FP+FN)
  c(TPR=TPR,FPR=FPR,PPV=PPV,F1=F1)
}

# =========================================================
# AGGREGATE NODE RESULTS
# =========================================================

num_reps <- length(unique(node$replicate))
num_locs <- 30
total_nodes <- num_locs * num_reps

node_summary <- node %>%
  group_by(method) %>%
  summarise(
    TP=sum(TP,na.rm=TRUE),
    FP=sum(FP,na.rm=TRUE),
    FN=sum(FN,na.rm=TRUE),
    .groups="drop"
  )

node_summary <- node_summary %>%
  mutate(TN = total_nodes - TP - FP - FN)

node_summary <- node_summary %>%
  rowwise() %>%
  mutate(vals=list(metrics(TP,FP,FN,TN))) %>%
  unnest_wider(vals)

# =========================================================
# AGGREGATE DOMAIN RESULTS
# =========================================================

t_length <- 100
total_domains <- num_locs * t_length * num_reps

domain_summary <- domain %>%
  group_by(method) %>%
  summarise(
    TP=sum(TP,na.rm=TRUE),
    FP=sum(FP,na.rm=TRUE),
    FN=sum(FN,na.rm=TRUE),
    .groups="drop"
  )

domain_summary <- domain_summary %>%
  mutate(TN = total_domains - TP - FP - FN)

domain_summary <- domain_summary %>%
  rowwise() %>%
  mutate(vals=list(metrics(TP,FP,FN,TN))) %>%
  unnest_wider(vals)

# =========================================================
# TYPE I ERROR AND POWER
# =========================================================

type1 <- node %>%
  filter(real_process==0) %>%
  summarise(rate = mean(tau < 0.05,na.rm=TRUE)) %>%
  pull(rate)

power <- node %>%
  filter(real_process==1) %>%
  summarise(rate = mean(tau < 0.05,na.rm=TRUE)) %>%
  pull(rate)

# =========================================================
# METHODS USED IN COMPARISON TABLE
# =========================================================

methods <- c(
  "DomainSelection_FocusGP_ISGM_TRUNC", # BMIM
  "FusedLasso",
  "Lasso",
  "FPCA_GP",
  "FLM_BS13",
  "FLM_BS23"
)

# =========================================================
# BEST VALUES FOR BOLDING
# =========================================================

best_TPR_d <- max(domain_summary$TPR,na.rm=TRUE)
best_FPR_d <- min(domain_summary$FPR,na.rm=TRUE)
best_PPV_d <- max(domain_summary$PPV,na.rm=TRUE)
best_F1_d  <- max(domain_summary$F1,na.rm=TRUE)

best_TPR_n <- max(node_summary$TPR,na.rm=TRUE)
best_FPR_n <- min(node_summary$FPR,na.rm=TRUE)
best_PPV_n <- max(node_summary$PPV,na.rm=TRUE)
best_F1_n  <- max(node_summary$F1,na.rm=TRUE)

# =========================================================
# TABLE 1: DOMAIN VS NO DOMAIN SELECTION
# =========================================================

bmim_domain <- node_summary[node_summary$method=="DomainSelection_FocusGP_ISGM_TRUNC",]
bmim_nodomain <- node_summary[node_summary$method=="DomainSelection_FocusGP_ISGM_FULL",]

cat("\n===== TABLE 1 ROWS =====\n\n")

cat("$\\mathbf{TPR}$ &",
    fmt(bmim_domain$TPR),"&",
    fmt(bmim_nodomain$TPR),"\\\\\n")

cat("$\\mathbf{FPR}$ &",
    fmt(bmim_domain$FPR),"&",
    fmt(bmim_nodomain$FPR),"\\\\\n")

cat("$\\mathbf{PPV}$ &",
    fmt(bmim_domain$PPV),"&",
    fmt(bmim_nodomain$PPV),"\\\\\n")

cat("$\\mathbf{F_{1}}$ &",
    fmt(bmim_domain$F1),"&",
    fmt(bmim_nodomain$F1),"\\\\\n")

# =========================================================
# TABLE 2: FULL METHOD COMPARISON
# =========================================================

cat("\n===== TABLE 2 ROWS =====\n\n")

# ------------------
# Domains
# ------------------

cat("\\parbox[t]{2mm}{\\multirow{4}{*}{\\rotatebox[origin=c]{90}{Domains}}}& $\\mathbf{TPR}$ &")

for(m in methods){
  val <- safe_val(domain_summary,m,"TPR")
  cat(" ",bold_max(val,best_TPR_d)," &")
}
cat("\\\\\n")

cat("&$\\mathbf{FPR}$ &")

for(m in methods){
  val <- safe_val(domain_summary,m,"FPR")
  cat(" ",bold_min(val,best_FPR_d)," &")
}
cat("\\\\\n")

cat("&$\\mathbf{PPV}$ &")

for(m in methods){
  val <- safe_val(domain_summary,m,"PPV")
  cat(" ",bold_max(val,best_PPV_d)," &")
}
cat("\\\\\n")

cat("&$\\mathbf{F_{1}}$ &")

for(m in methods){
  val <- safe_val(domain_summary,m,"F1")
  cat(" ",bold_max(val,best_F1_d)," &")
}
cat("\\\\\n")

# ------------------
# Nodes
# ------------------

cat("\\Xhline{\\arrayrulewidth}\n")

cat("\\parbox[t]{2mm}{\\multirow{4}{*}{\\rotatebox[origin=c]{90}{Nodes}}}& $\\mathbf{TPR}$ &")

for(m in methods){
  val <- safe_val(node_summary,m,"TPR")
  cat(" ",bold_max(val,best_TPR_n)," &")
}
cat("\\\\\n")

cat("&$\\mathbf{FPR}$ &")

for(m in methods){
  val <- safe_val(node_summary,m,"FPR")
  cat(" ",bold_min(val,best_FPR_n)," &")
}
cat("\\\\\n")

cat("&$\\mathbf{PPV}$ &")

for(m in methods){
  val <- safe_val(node_summary,m,"PPV")
  cat(" ",bold_max(val,best_PPV_n)," &")
}
cat("\\\\\n")

cat("&$\\mathbf{F_{1}}$ &")

for(m in methods){
  val <- safe_val(node_summary,m,"F1")
  cat(" ",bold_max(val,best_F1_n)," &")
}
cat("\\\\\n")

# ------------------
# Overall
# ------------------

cat("\\Xhline{\\arrayrulewidth}\n")

cat("\\parbox[t]{2mm}{\\multirow{2}{*}{\\rotatebox[origin=c]{90}{Overall}}}&$P(\\text{Reject }H_{0}|H_{0})$ &",
    fmt(type1)," & - & - & - & - & -\\\\\n")

cat("&$P(\\text{Reject }H_{0}|H_{1})$ &",
    fmt(power)," & - & - & - & - & -\\\\\n")