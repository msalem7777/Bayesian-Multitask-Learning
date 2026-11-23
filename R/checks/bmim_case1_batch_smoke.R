# ---------------------------------------------------------------------------
# BMIM case-1 batch CLI smoke test: no simulations and no output files.
# Run from the R/ directory:
#   Rscript --vanilla -e "source('checks/bmim_case1_batch_smoke.R')"
# ---------------------------------------------------------------------------
script <- "bmim_case1_batch.R"
if (!file.exists(script)) {
  stop("Run the batch smoke from the repository's R/ directory.")
}
result <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", script, "--self-test"),
                  stdout = TRUE, stderr = TRUE)
status <- attr(result, "status")
if (!is.null(status) && status != 0L) {
  cat(paste(result, collapse = "\n"), "\n")
  stop("BMIM batch CLI self-test failed: status ", status)
}
if (!any(grepl("BMIM CASE 1 BATCH RUNNER SELF-TEST: PASS", result,
               fixed = TRUE))) {
  stop("Expected success message missing from batch CLI self-test.")
}
cat("BMIM CASE 1 BATCH SMOKE: PASS\n")
