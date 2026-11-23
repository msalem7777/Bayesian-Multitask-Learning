# =========================================================
# generate_latex_tables.R
#
# Reads c2_domain_results_10.csv and c2_node_results_10.csv,
# computes TPR / FPR / PPV / F1 using the same logic as
# evaluation_metrics.R, and outputs two LaTeX tables:
#
#  Table 1 – Domain Selection vs No Domain Selection (node-level)
#             ISGM_FULL  vs  ISGM_TRUNC
#
#  Table 2 – Full comparison: BMIM (ISGM_TRUNC) vs all baselines
#             Domain rows + Node rows + Overall significance rows
# =========================================================

source("C:/My Files/repos/Bayesian-Multitask-Learning/R/evaluation_metrics.R")   # provides selection_metrics(), etc.

domain_df <- read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_domain_results_10_LN_newD.csv", stringsAsFactors = FALSE)
node_df   <- read.csv("C:/My Files/repos/Bayesian-Multitask-Learning/R/c2_node_results_10_LN_newD.csv",   stringsAsFactors = FALSE)

# ── num_locs inferred from data ────────────────────────────────────────────
num_locs <- max(domain_df$location, na.rm = TRUE)   # e.g. 10


# =========================================================
# HELPER: aggregate TP / FP / FN across replicates then
#         compute macro TPR, FPR, PPV, F1
#
# For FPR we need TN.  In the node table the universe is
# num_locs per replicate; in the domain table the universe
# is t_length time-points per location per replicate.
# Both tables already store TP / FP / FN directly.
# TN = (universe_size) - TP - FP - FN
# =========================================================

# ── NODE-LEVEL metrics ─────────────────────────────────────────────────────
# Each row in node_df is already one replicate (TP/FP/FN summed across num_locs).
# Compute TPR/FPR/PPV/F1 per replicate then average across replicates.
node_metrics <- function(df, method_name, scenario) {
  sub <- df[df$method == method_name & df$real_process == scenario, ]
  
  rates <- lapply(seq_len(nrow(sub)), function(i) {
    TP <- sub$TP[i];  FP <- sub$FP[i];  FN <- sub$FN[i]
    TN <- num_locs - TP - FP - FN
    TPR <- if (!is.na(TP + FN) && (TP + FN) > 0) TP / (TP + FN) else NA
    FPR <- if (!is.na(FP + TN) && (FP + TN) > 0) FP / (FP + TN) else NA
    PPV <- if (!is.na(TP + FP) && (TP + FP) > 0) TP / (TP + FP) else NA
    F1  <- if (!is.na(PPV + TPR) && (PPV + TPR) > 0) 2*PPV*TPR/(PPV+TPR) else NA
    c(TPR = TPR, FPR = FPR, PPV = PPV, F1 = F1)
  })
  
  rates_mat <- do.call(rbind, rates)
  list(
    TPR = mean(rates_mat[, "TPR"], na.rm = TRUE),
    FPR = mean(rates_mat[, "FPR"], na.rm = TRUE),
    PPV = mean(rates_mat[, "PPV"], na.rm = TRUE),
    F1  = mean(rates_mat[, "F1"],  na.rm = TRUE)
  )
}

# ── DOMAIN-LEVEL metrics ────────────────────────────────────────────────────
# Each row in domain_df is one node (location) within one replicate.
# Step 1: sum TP/FP/FN across all P nodes within each replicate.
# Step 2: compute TPR/FPR/PPV/F1 per replicate (need TN: t_length*P - TP-FP-FN).
# Step 3: average the per-replicate rates across replicates.
# This mirrors how node_metrics() works but at the time-point level.
t_length <- 100   # matches cfg$t_length = 100 in simulation_runs_new_v2.R

domain_metrics <- function(df, method_name, scenario) {
  sub <- df[df$method == method_name & df$real_process == scenario, ]
  
  # Sum TP/FP/FN per replicate across all nodes, ignoring NA rows
  rep_ids <- unique(sub$replicate)
  
  rates <- lapply(rep_ids, function(rid) {
    rows <- sub[sub$replicate == rid, ]
    TP <- sum(rows$TP, na.rm = TRUE)
    FP <- sum(rows$FP, na.rm = TRUE)
    FN <- sum(rows$FN, na.rm = TRUE)
    TN <- num_locs * t_length - TP - FP - FN
    TPR <- if ((TP + FN) > 0) TP / (TP + FN) else NA
    FPR <- if ((FP + TN) > 0) FP / (FP + TN) else NA
    PPV <- if ((TP + FP) > 0) TP / (TP + FP) else NA
    F1  <- if (!is.na(PPV + TPR) && (PPV + TPR) > 0) 2*PPV*TPR/(PPV+TPR) else NA
    c(TPR = TPR, FPR = FPR, PPV = PPV, F1 = F1)
  })
  
  rates_mat <- do.call(rbind, rates)
  list(
    TPR = mean(rates_mat[, "TPR"], na.rm = TRUE),
    FPR = mean(rates_mat[, "FPR"], na.rm = TRUE),
    PPV = mean(rates_mat[, "PPV"], na.rm = TRUE),
    F1  = mean(rates_mat[, "F1"],  na.rm = TRUE)
  )
}

# ── MEAN IoU ──────────────────────────────────────────────────────────────
# Average per-replicate mean IoU (first average across nodes ignoring NA,
# then average those replicate means across replicates)
mean_iou <- function(df, method_name, scenario) {
  sub <- df[df$method == method_name & df$real_process == scenario, ]
  per_rep <- tapply(sub$IoU, sub$replicate, function(x) mean(x, na.rm = TRUE))
  mean(per_rep, na.rm = TRUE)
}

# ── OVERALL significance ───────────────────────────────────────────────────
# P(Reject H0 | H0): proportion of S0 replicates where tau < 0.05
# P(Reject H0 | H1): proportion of S1 replicates where tau < 0.05
overall_metrics <- function(df, method_name) {
  s0 <- df[df$method == method_name & df$real_process == 0, ]
  s1 <- df[df$method == method_name & df$real_process == 1, ]
  
  type1 <- if (nrow(s0) > 0) mean(s0$tau < 0.05, na.rm = TRUE) else NA
  power <- if (nrow(s1) > 0) mean(s1$tau < 0.05, na.rm = TRUE) else NA
  
  list(type1 = type1, power = power)
}

# ── Formatting helpers ─────────────────────────────────────────────────────
fmt <- function(x, digits = 3) {
  if (is.na(x)) return("-")
  formatC(round(x, digits), format = "f", digits = digits)
}

# Bold the better value between two (lower is better for FPR, higher for rest)
bold_better <- function(a, b, lower_is_better = FALSE) {
  if (is.na(a) || is.na(b)) return(c(fmt(a), fmt(b)))
  if (lower_is_better) {
    if (a <= b) return(c(paste0("\\bf ", fmt(a)), fmt(b)))
    else        return(c(fmt(a), paste0("\\bf ", fmt(b))))
  } else {
    if (a >= b) return(c(paste0("\\bf ", fmt(a)), fmt(b)))
    else        return(c(fmt(a), paste0("\\bf ", fmt(b))))
  }
}

# Bold the maximum value among a vector (lower for FPR)
bold_max_vec <- function(vals, lower_is_better = FALSE) {
  out <- sapply(vals, fmt)
  valid <- which(!is.na(vals))
  if (length(valid) == 0) return(out)
  if (lower_is_better) best <- valid[which.min(vals[valid])]
  else                  best <- valid[which.max(vals[valid])]
  out[best] <- paste0("\\bf ", fmt(vals[best]))
  out
}


# =========================================================
# TABLE 1 – Domain Selection vs No Domain Selection (Nodes)
#
# Column 1 (Domain Selection):    A3_Domain_ISGM (node_df)
#                                  = DomainSelection_FocusGP_ISGM_TRUNC (domain_df)
# Column 2 (No Domain Selection):  A4_ISGM_Full   (node_df)
#                                  = DomainSelection_FocusGP_ISGM_FULL  (domain_df)
# Both evaluated on node (location) level under S1.
# =========================================================

m_full  <- node_metrics(node_df, "A3_Domain_ISGM", scenario = 1)  # Domain Selection
m_trunc <- node_metrics(node_df, "A4_ISGM_Full",  scenario = 1)  # No Domain Selection

tpr_row <- bold_better(m_full$TPR,  m_trunc$TPR,  lower_is_better = FALSE)
fpr_row <- bold_better(m_full$FPR,  m_trunc$FPR,  lower_is_better = TRUE)
ppv_row <- bold_better(m_full$PPV,  m_trunc$PPV,  lower_is_better = FALSE)
f1_row  <- bold_better(m_full$F1,   m_trunc$F1,   lower_is_better = FALSE)

table1 <- paste0(
  "\\begin{table}[H]\n",
  "\\centering\n",
  "\\begin{tabular}{|c|c|c|c|}\n",
  "\\hline\n",
  "    &&Domain Selection &  No Domain Selection\\\\\n",
  "    \\hline\n",
  "    \\parbox[t]{2mm}{\\multirow{4}{*}{\\rotatebox[origin=c]{90}{Nodes}}}",
  "&  $\\mathbf{TPR}$ & ", tpr_row[1], " & ", tpr_row[2], " \\\\\n",
  "    &$\\mathbf{FPR}$ & ", fpr_row[1], " & ", fpr_row[2], " \\\\\n",
  "    &$\\mathbf{PPV}$ & ", ppv_row[1], " & ", ppv_row[2], " \\\\\n",
  "    &$\\mathbf{F_{1}}$ & ", f1_row[1], " & ", f1_row[2], " \\\\\n",
  "    \\hline \n",
  "\\end{tabular}\n",
  "\\caption{Impact of BMIM domain selection on task of node selection when N=200 and P=10}\n",
  "\\label{C2_BMIM_P10_n200}\n",
  "\\end{table}"
)

cat("\n\n% ── TABLE 1 ─────────────────────────────────────────\n")
cat(table1)


# =========================================================
# TABLE 2 – Full comparison
#
# Columns: BMIM (=ISGM_TRUNC / A3_Domain_ISGM) | FusedLasso | Lasso |
#          FPCA-GP | FLM-BS(13) | FLM-BS(23)
#
# Domain rows   → domain_df,  S1 only
# Node rows     → node_df,    S1 only
# Overall rows  → node_df tau, S0 + S1
#
# Method mapping:
#   BMIM domain  → DomainSelection_FocusGP_ISGM_TRUNC  (domain_df)
#   BMIM node    → A3_Domain_ISGM                       (node_df)
#   FusedLasso   → FusedLasso
#   Lasso        → Lasso
#   FPCA-GP      → GP_FPCA / FPCA_GP
#   FLM-BS(13)   → BSpline_13 / FLM_BS13
#   FLM-BS(23)   → BSpline_23 / FLM_BS23
# =========================================================

# ── Domain rows ────────────────────────────────────────────────────────────
dm_BMIM  <- domain_metrics(domain_df, "DomainSelection_FocusGP_ISGM_TRUNC", 1)  # A3 in domain CSV
dm_FL    <- domain_metrics(domain_df, "FusedLasso",  1)
dm_L     <- domain_metrics(domain_df, "Lasso",       1)
dm_FPCA  <- list(TPR=NA, FPR=NA, PPV=NA, F1=NA)   # FPCA-GP has no domain sel
dm_BS13  <- domain_metrics(domain_df, "FLM_BS13",    1)
dm_BS23  <- domain_metrics(domain_df, "FLM_BS23",    1)

iou_BMIM <- mean_iou(domain_df, "DomainSelection_FocusGP_ISGM_TRUNC", 1)
iou_FL   <- mean_iou(domain_df, "FusedLasso", 1)
iou_L    <- mean_iou(domain_df, "Lasso",      1)
iou_BS13 <- mean_iou(domain_df, "FLM_BS13",   1)
iou_BS23 <- mean_iou(domain_df, "FLM_BS23",   1)

dom_TPR <- bold_max_vec(c(dm_BMIM$TPR, dm_FL$TPR, dm_L$TPR, NA, dm_BS13$TPR, dm_BS23$TPR), FALSE)
dom_FPR <- bold_max_vec(c(dm_BMIM$FPR, dm_FL$FPR, dm_L$FPR, NA, dm_BS13$FPR, dm_BS23$FPR), TRUE)
dom_PPV <- bold_max_vec(c(dm_BMIM$PPV, dm_FL$PPV, dm_L$PPV, NA, dm_BS13$PPV, dm_BS23$PPV), FALSE)
dom_F1  <- bold_max_vec(c(dm_BMIM$F1,  dm_FL$F1,  dm_L$F1,  NA, dm_BS13$F1,  dm_BS23$F1),  FALSE)
dom_IoU <- bold_max_vec(c(iou_BMIM,    iou_FL,    iou_L,    NA, iou_BS13,    iou_BS23),     FALSE)

# ── Node rows ──────────────────────────────────────────────────────────────
nm_BMIM  <- node_metrics(node_df, "A3_Domain_ISGM", 1)
nm_FL    <- node_metrics(node_df, "FusedLasso",      1)
nm_L     <- node_metrics(node_df, "Lasso",           1)
nm_FPCA  <- node_metrics(node_df, "GP_FPCA",         1)
nm_BS13  <- node_metrics(node_df, "BSpline_13",      1)
nm_BS23  <- node_metrics(node_df, "BSpline_23",      1)

nod_TPR <- bold_max_vec(c(nm_BMIM$TPR, nm_FL$TPR, nm_L$TPR, nm_FPCA$TPR, nm_BS13$TPR, nm_BS23$TPR), FALSE)
nod_FPR <- bold_max_vec(c(nm_BMIM$FPR, nm_FL$FPR, nm_L$FPR, nm_FPCA$FPR, nm_BS13$FPR, nm_BS23$FPR), TRUE)
nod_PPV <- bold_max_vec(c(nm_BMIM$PPV, nm_FL$PPV, nm_L$PPV, nm_FPCA$PPV, nm_BS13$PPV, nm_BS23$PPV), FALSE)
nod_F1  <- bold_max_vec(c(nm_BMIM$F1,  nm_FL$F1,  nm_L$F1,  nm_FPCA$F1,  nm_BS13$F1,  nm_BS23$F1),  FALSE)

# ── Overall significance (only BMIM has tau) ───────────────────────────────
ov_BMIM <- overall_metrics(node_df, "A3_Domain_ISGM")
type1_str <- fmt(ov_BMIM$type1)
power_str <- fmt(ov_BMIM$power)

table2 <- paste0(
  "\\begin{table}[H]\n",
  "    \\centering\n",
  "\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n",
  "    \\hline\\hline\n",
  "        && BMIM & Fused Lasso & Lasso & FPCA-GP & FLM-BS(13) & FLM-BS(23)\\\\\n",
  "    \\hline\n",
  # Domain rows
  "    \\parbox[t]{2mm}{\\multirow{5}{*}{\\rotatebox[origin=c]{90}{Domains}}}",
  "& $\\mathbf{TPR}$ & ", dom_TPR[1], " & ", dom_TPR[2], " & ", dom_TPR[3],
  " & ", dom_TPR[4], "& ", dom_TPR[5], " & ", dom_TPR[6], "\\\\\n",
  "    &$\\mathbf{FPR}$ & ", dom_FPR[1], " & ", dom_FPR[2], "  & ", dom_FPR[3],
  " & ", dom_FPR[4], " & ", dom_FPR[5], " & ", dom_FPR[6], "\\\\\n",
  "    &$\\mathbf{PPV}$ & ", dom_PPV[1], " & ", dom_PPV[2], " & ", dom_PPV[3],
  " & ", dom_PPV[4], " & ", dom_PPV[5], " & ", dom_PPV[6], "\\\\\n",
  "    &$\\mathbf{F_{1}}$ &  ", dom_F1[1], " & ", dom_F1[2], " & ", dom_F1[3],
  " & ", dom_F1[4], " & ", dom_F1[5], " & ", dom_F1[6], "\\\\\n",
  "    &$\\mathbf{IoU}$ & ", dom_IoU[1], " & ", dom_IoU[2], " & ", dom_IoU[3],
  " & ", dom_IoU[4], " & ", dom_IoU[5], " & ", dom_IoU[6], "\\\\\n",
  "    \\Xhline{\\arrayrulewidth}\n",
  "    \\hline\n",
  # Node rows
  "    \\parbox[t]{2mm}{\\multirow{4}{*}{\\rotatebox[origin=c]{90}{Nodes}}}",
  "& $\\mathbf{TPR}$ & ", nod_TPR[1], " & ", nod_TPR[2], " & ", nod_TPR[3],
  " & ", nod_TPR[4], " & ", nod_TPR[5], " & ", nod_TPR[6], "\\\\\n",
  "    &$\\mathbf{FPR}$ &  ", nod_FPR[1], " & ", nod_FPR[2], "  & ", nod_FPR[3],
  " & ", nod_FPR[4], " & ", nod_FPR[5], " & ", nod_FPR[6], "\\\\\n",
  "    &$\\mathbf{PPV}$ & ", nod_PPV[1], " & ", nod_PPV[2], " & ", nod_PPV[3],
  " & ", nod_PPV[4], "& ", nod_PPV[5], " & ", nod_PPV[6], "\\\\\n",
  "    &$\\mathbf{F_{1}}$ & ", nod_F1[1], " & ", nod_F1[2], "  & ", nod_F1[3],
  " & ", nod_F1[4], "& ", nod_F1[5], " & ", nod_F1[6], "\\\\\n",
  "    \\Xhline{\\arrayrulewidth}\n",
  # Overall rows
  "    \\parbox[t]{2mm}{\\multirow{2}{*}{\\rotatebox[origin=c]{90}{Overall}}}",
  "&$P(\\text{Reject }H_{0}|H_{0})$ & ", type1_str,
  " & - & - & - & - & -\\\\\n",
  "    &$P(\\text{Reject }H_{0}|H_{1})$ & ", power_str,
  " & - & - & - & - & -\\\\\n",
  "    \\hline \n",
  "\\end{tabular}\n",
  "\\caption{Comparison against other valid methodologies under Case 2: BMIM, ",
  "Fused Lasso, Lasso, FPCA-GP, FLM-BS(13): FLM-BSpline (B=13), FLM-BS(23): FLM-BSpline (B=23). ",
  "The table shows the chosen models' (1) domain selection performance metrics, ",
  "(2) node selection performance metrics, (3) overall significance under the null (Type-I error), ",
  "and (4) overall significance under the alternative (Power) when N=200 and P=10}\n",
  "\\label{C2_All_P10_n200}\n",
  "\\end{table}"
)

cat("\n\n% ── TABLE 2 ─────────────────────────────────────────\n")
cat(table2)
cat("\n")