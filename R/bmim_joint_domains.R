# ============================================================================
# bmim_joint_domains.R -- modular loader for the JOINT DOMAIN EXTENSION.
# Legacy fixed-domain defaults remain unchanged; scalar-aware callbacks
# are available in bmim_gp_updates.R.
# First source bmim_gp_updates.R and bmim_global_gp_test.R (as needed).
# Then source THIS loader; it installs all standalone joint-domain modules.
# ============================================================================
.bmim_joint_files <- Filter(
  Negate(is.null),
  lapply(sys.frames(), function(frame) frame$ofile)
)
if (length(.bmim_joint_files) == 0L) {
  stop("Source bmim_joint_domains.R from its actual file path.")
}
.bmim_joint_root <- dirname(normalizePath(tail(.bmim_joint_files,1L)[[1L]]))
for (.bmim_joint_name in c("bmim_domain_distances.R",
                           "bmim_scalar_effects.R",
                           "bmim_vb_prior.R",
                           "bmim_domain_updates.R",
                           "bmim_joint_sampler.R")) {
  source(file.path(.bmim_joint_root,.bmim_joint_name), local=TRUE)
}
rm(.bmim_joint_files,.bmim_joint_root,.bmim_joint_name)
