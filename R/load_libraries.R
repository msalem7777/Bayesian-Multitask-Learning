# =========================================================
# load_libraries.R
# Shared library loader — source this at the top of any
# script in the project before sourcing other modules.
# =========================================================

library(data.table)
library(invgamma)
library(hetGP)
library(mlegp)
library(DiceKriging)
library(infotheo)   # mutinformation(), discretize()
library(mvtnorm)    # rmvnorm()
library(laGP)       # newGPsep(), mleGPsep(), predGPsep(), deleteGPsep(), darg()
library(plgp)       # covar.sep(), distance()
library(lhs)        # randomLHS()
library(tgp)        # bgp(), tgp.default.params()
library(MASS)
library(Matrix)
library(astsa)
library(mgcv)
library(fda)        # create.bspline.basis(), smooth.basisPar()
library(fdapace)    # FPCA()
library(extlasso)   # fusedlasso()
library(glmnet)     # cv.glmnet()
library(refund)
library(ClusterR)
library(ggplot2)