# BMIM validation entry point: execute from the repository R directory.
# Smoke checks are quick; the seven-configuration numerical integration and
# four-chain exploration audit is deliberately more demanding.
# PowerShell: & $rs -e "source('checks/run_all_tests.R')"
files <- c('bmim_gp_smoke.R',
           'bmim_global_gp_smoke.R',
           'bmim_global_gp_eight_state_toy.R',
           'bmim_toy_exact_validation.R')
source_files <- Filter(Negate(is.null),
                       lapply(sys.frames(), function(frame) frame$ofile))
if (!length(source_files)) stop('Source run_all_tests.R from its file.')
check_dir <- dirname(normalizePath(tail(source_files, 1L)[[1L]]))
cat('BMIM TEST SUITE\n')
for (filename in files) {
  path <- file.path(check_dir, filename)
  if (!file.exists(path)) stop('Missing test: ', path)
  cat('\n=== ', filename, ' ===\n', sep='')
  source(path, local=new.env(parent=globalenv()))
}
cat('\nBMIM ALL FOUR VALIDATION SCRIPTS: PASS\n')
