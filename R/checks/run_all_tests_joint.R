# End-to-end BMIM validation entry point. Run from R/:
#   & $rs -e "source('checks/run_all_tests_joint.R')"
# Existing validated checks are kept intact; the two new joint-domain
# suites are additive and cannot affect fixed-domain behavior.
files <- Filter(Negate(is.null),
                lapply(sys.frames(),function(frame) frame$ofile))
if (!length(files)) stop("Source the suite from its actual .R file.")
check_dir <- dirname(normalizePath(tail(files,1L)[[1L]]))
cat("=== LEGACY VALIDATED FOUR-CHECK SUITE ===\n")
source(file.path(check_dir,"run_all_tests.R"), local=new.env(parent=globalenv()))
for (name in c("bmim_joint_domain_smoke.R",
               "bmim_joint_domain_gp_smoke.R")) {
  cat("\n=== ",name," ===\n",sep="")
  source(file.path(check_dir,name),local=new.env(parent=globalenv()))
}
cat("\nBMIM ALL SIX VALIDATION SCRIPTS: PASS\n")
