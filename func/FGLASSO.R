###############################################################
#                        FGLASSO
# This script implements the Functional Graphical Lasso (FGLASSO)
# algorithm using Functional Principal Component Analysis (FPCA).
# It estimates the sparse structure of covariance among functional data.
###############################################################

library(fda)

# ------------------------------------------------------------
# Function: FGLASSO
# Purpose: Perform Functional Graphical Lasso estimation.
# Inputs:
#   h       : A 3D array of functional data with dimensions:
#             n (number of samples) x tau (number of time points) x p (number of variables)
#   M       : Number of B-spline basis functions.
#   L       : Number of lambda values to generate along the regularization path.
#   obs.time: Vector of observation time points.
# Output:
#   A list of support matrices (each p x p) corresponding to each lambda value.
# Process:
#   - Convert each variable's functional data into a matrix using a B-spline basis.
#   - Perform FPCA to obtain principal component scores.
#   - Center the FPCA scores and compute their covariance matrix.
#   - Compute the maximum lambda value using lambda.sup.global.gY.
#   - Generate a sequence of lambda values.
#   - For each lambda, compute the support matrix using the ProxAlg_FGM function.
# ------------------------------------------------------------
FGLASSO <- function(h, M, L, obs.time) {
  p <- dim(h)[3]
  tau <- dim(h)[2]
  n <- dim(h)[1]
  fpc.score <- numeric(0)
  
  for (j in 1:p) {
    obs.val.matrix <- matrix(0, nrow = tau, ncol = n)
    for (i in 1:n) {
      obs.val.vec <- as.vector(h[i, , j])
      obs.val.matrix[, i] <- obs.val.vec
    }
    bspline.basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = M)
    fd.object.array <- Data2fd(argvals = obs.time, y = obs.val.matrix, basisobj = bspline.basis)
    # Perform FPCA process
    fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm = M)$scores)
  }
  
  fpc.score.cen <- scale(fpc.score, center = TRUE, scale = FALSE)
  est.cov.pc <- (t(fpc.score.cen) %*% fpc.score.cen) / (n - 1)
  lambda.max.gY <- lambda.sup.global.gY(h1, M)  # Note: 'h1' should be defined or passed appropriately
  lambdas.gY <- seq(lambda.max.gY, 0, length.out = L)
  
  G.fglasso <- list()
  for (l in 1:L) {
    # Compute the p x p support matrix for the l-th lambda value
    G.fglasso[[l]] <- ProxAlg_FGM(est.cov.pc, p, M, lambdas.gY[l])$Support
  }
  return(G.fglasso)
}
