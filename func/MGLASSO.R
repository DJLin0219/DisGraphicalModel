###############################################################
#                       MGLASSO
# This script implements the Functional Graphical Lasso (MGLASSO)
# algorithm for estimating sparse graphical structures in functional data.
# It includes functions to compute kernel representations, covariance
# and correlation matrices, and various implementations of the MGLASSO algorithm.
###############################################################

# ------------------------------------------------------------
# Function: RedKer2
# Purpose: Compute a reduced kernel representation by performing eigenvalue decomposition
#          on matrix G and selecting the top C eigencomponents.
# Inputs:
#   G : A square matrix.
#   C : Number of top eigencomponents to retain.
# Output:
#   A list with:
#     SCORE: Matrix of scores computed as the product of the top eigenvectors and the square roots of their eigenvalues.
#     d    : The largest eigenvalue.
# ------------------------------------------------------------
RedKer2 <- function(G, C) {
  # Perform eigenvalue decomposition
  eig <- eigen(G)
  U <- eig$vectors
  D <- eig$values
  
  # Sort eigenvalues in descending order
  sorted_indices <- order(D, decreasing = TRUE)
  D_sorted <- D[sorted_indices]
  U_sorted <- U[, sorted_indices]
  
  # Select top C eigenvectors and compute scores
  COEF <- U_sorted[, 1:C]
  SCORE <- COEF %*% diag(sqrt(D_sorted[1:C]))
  
  # Extract the largest eigenvalue
  d <- D_sorted[1]
  
  list(SCORE = SCORE, d = d)
}

# ------------------------------------------------------------
# Function: mppower
# Purpose: Compute the Moore-Penrose power of a matrix using its eigen decomposition.
# Inputs:
#   matrix : A square matrix.
#   power  : The exponent to which the eigenvalues are raised.
#   ignore : Threshold to ignore small eigenvalues.
# Output:
#   The matrix raised to the specified power.
# ------------------------------------------------------------
mppower = function(matrix, power, ignore) {
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval) > ignore])
  if (m == 1) {
    tmp = (eval[1:m]^power) * evec[, 1:m] %*% t(evec[, 1:m])
  } else {
    tmp = evec[, 1:m] %*% diag(eval[1:m]^power) %*% t(evec[, 1:m])
  }
  return(tmp)
}

# ------------------------------------------------------------
# Function: cov_to_cor
# Purpose: Convert a covariance matrix representation into a correlation-like form using ridge regularization.
# Inputs:
#   A     : A matrix (typically a score matrix).
#   ridge : Regularization parameter.
#   n     : Sample size.
# Output:
#   A transformed matrix computed using the inverse square root of a regularized version of A'A.
# ------------------------------------------------------------
cov_to_cor <- function(A, ridge, n) {
  ni <- ncol(A)
  ATA <- t(A) %*% A
  ATAridge <- ATA + n * ridge * qmat(ni)
  # The following line computes a regularized version of ATA for inversion
  ATAinv <- t(A) %*% A + n * ridge * eigen(ATA)$values[1] * qmat(ni)
  return(t(mppower(ATAinv, power = -1/2, ignore = 1e-6) %*% t(A)))
}

# ------------------------------------------------------------
# Function: Cov.new
# Purpose: Compute a new covariance matrix from functional data.
# Inputs:
#   N.t   : A vector of counts for each variable.
#   xvec  : A 3D array of functional data (n x tau x p).
#   ridge : Regularization parameter.
#   M     : Number of components to retain (default is 5).
# Output:
#   A symmetric covariance matrix.
# Process:
#   - Compute a kernel matrix for each variable using gramx.new and nzmat functions.
#   - Apply RedKer2 to obtain the score representation.
#   - Assemble the overall covariance matrix from the concatenated scores.
# ------------------------------------------------------------
Cov.new <- function(N.t, xvec, ridge, M = 5) {
  n = dim(xvec)[1] / max(N.t)
  p = dim(xvec)[3]
  eps = 10^(-7)
  ns = dim(xvec)[2]
  kt = diag(simpson(ns))
  # Initialize KX matrix to store scores
  KX <- matrix(0, n, p * M)
  for (j in 1:p) {
    Gj <- cgram(gramx.new(N.t[j], nzmat(xvec[, , j]), kt))
    KX[,(1 + (j - 1) * M):(j * M)] <- RedKer2(Gj, C = M)$SCORE / sqrt(n)
  }
  Cov <- t(KX) %*% KX
  return((Cov + t(Cov)) / 2)
}

# ------------------------------------------------------------
# Function: Corr.new
# Purpose: Compute a new correlation matrix from functional data.
# Inputs:
#   N.t   : A vector of counts for each variable.
#   xvec  : A 3D array of functional data.
#   ridge : Regularization parameter.
#   M     : Number of components to retain (default is 5).
# Output:
#   A symmetric correlation matrix.
# Process:
#   - For each variable, compute a kernel representation using gramx.new and nzmat.
#   - Transform the score matrix to a correlation-like form using cov_to_cor.
#   - Assemble the overall correlation matrix.
# ------------------------------------------------------------
Corr.new <- function(N.t, xvec, ridge, M = 5) {
  n = dim(xvec)[1] / max(N.t)
  p = dim(xvec)[3]
  eps = 10^(-7)
  ns = dim(xvec)[2]
  kt = diag(simpson(ns))
  # Initialize KX matrix to store transformed scores
  KX <- matrix(0, n, p * M)
  for (j in 1:p) {
    Gj <- cgram(gramx.new(N.t[j], nzmat(xvec[, , j]), kt))
    KX[,(1 + (j - 1) * M):(j * M)] <- cov_to_cor(RedKer2(Gj, C = M)$SCORE, ridge, n)
  }
  Corr <- t(KX) %*% KX
  return((Corr + t(Corr)) / 2)
}

# ------------------------------------------------------------
# Function: Corr.new2
# Purpose: Compute a new correlation matrix using an alternative grouping method.
# Inputs:
#   xvec  : A 3D array of functional data.
#   ridge : Regularization parameter.
#   Ndim  : A vector specifying dimensions for grouping variables.
#   p     : Number of variables.
#   M     : Number of components to retain (default is 5).
# Output:
#   A symmetric correlation matrix.
# Process:
#   - Group the data according to Ndim.
#   - Compute the kernel representation using gramx.new2.
#   - Transform the score matrix to a correlation-like form using cov_to_cor.
#   - Assemble the overall correlation matrix.
# ------------------------------------------------------------
Corr.new2 <- function(xvec, ridge, Ndim, p, M = 5) {
  n = dim(xvec)[1]
  ns = dim(xvec)[2]
  eps = 10^(-7)
  kt = diag(simpson(ns))
  KX <- matrix(0, n, p * M)
  for (j in 1:p) {
    index <- (sum(Ndim[1:(j - 1)]) + 1):(sum(Ndim[1:j]))
    Gj <- cgram(gramx.new2(xvec[, , index], kt))
    KX[,(1 + (j - 1) * M):(j * M)] <- cov_to_cor(RedKer2(Gj, C = M)$SCORE, ridge, n)
  }
  Corr <- t(KX) %*% KX
  return((Corr + t(Corr)) / 2)
}

# ------------------------------------------------------------
# Function: Cov.new2
# Purpose: Compute a new covariance matrix using an alternative grouping method.
# Inputs:
#   xvec  : A 3D array of functional data.
#   ridge : Regularization parameter.
#   Ndim  : A vector specifying dimensions for grouping variables.
#   p     : Number of variables.
#   M     : Number of components to retain (default is 5).
# Output:
#   A symmetric covariance matrix.
# Process:
#   - Group the data according to Ndim.
#   - Compute the kernel representation using gramx.new2 and obtain scores via RedKer2.
#   - Assemble the overall covariance matrix.
# ------------------------------------------------------------
Cov.new2 <- function(xvec, ridge, Ndim, p, M = 5) {
  n = dim(xvec)[1]
  ns = dim(xvec)[2]
  eps = 10^(-7)
  kt = diag(simpson(ns))
  KX <- matrix(0, n, p * M)
  for (j in 1:p) {
    index <- (sum(Ndim[1:(j - 1)]) + 1):(sum(Ndim[1:j]))
    Gj <- cgram(gramx.new2(xvec[, , index], kt))
    KX[,(1 + (j - 1) * M):(j * M)] <- RedKer2(Gj, C = M)$SCORE / sqrt(n)
  }
  Cov <- t(KX) %*% KX
  return((Cov + t(Cov)) / 2)
}

# ------------------------------------------------------------
# Function: MGLASSO
# Purpose: Perform Functional Graphical Lasso estimation using the computed correlation matrix to estimate graph structure when the data
#          has dimensions (n*N) x tau x p.
# Inputs:
#   data  : A 3D array of functional data.
#   n     : Sample size.
#   N.t   : A vector of counts for each variable.
#   M     : Number of components (default is 5).
#   L     : Number of lambda values along the regularization path (default is 100).
#   ridge : Regularization parameter.
# Output:
#   A list of support matrices (each p x p) corresponding to each lambda value.
# Process:
#   - Compute the correlation matrix using Corr.new.
#   - Generate a sequence of lambda values from max(est.cov.pc) to 0.
#   - For each lambda, apply ProxAlg_FGM to compute the support matrix.
# ------------------------------------------------------------
MGLASSO <- function(data, n, N.t, M = 5, L = 100, ridge = ridge) {
  p <- dim(data)[3]
  tau <- dim(data)[2]
  est.cov.pc <- Corr.new(N.t = N.t, xvec = data, ridge = ridge, M = 5)
  lambda.cmcg.max <- max(est.cov.pc)
  lambda.cmcg <- seq(lambda.cmcg.max, 0, length.out = L)
  G.fglasso <- list()
  for (j in 1:L) {
    G.fglasso[[j]] <- ProxAlg_FGM(S = est.cov.pc, p = p, M = n, gamma = lambda.cmcg[j],
                                  Eta = "Auto", n.iteration = 10)$Support  # p x p support matrix
    print(j)
  }
  return(G.fglasso)
}

# ------------------------------------------------------------
# Function: MGLASSO2
# Purpose: Perform Functional Graphical Lasso estimation using the computed correlation matrix to estimate graph structure when the data
#          has dimensions (n) x tau x (N1 + ... + Np).
# Inputs:
#   data  : A 3D array of functional data.
#   p     : Number of variables.
#   Ndim  : A vector specifying dimensions for grouping variables.
#   M     : Number of components (default is 5).
#   L     : Number of lambda values along the regularization path (default is 100).
#   ridge : Regularization parameter.
# Output:
#   A list of support matrices (each p x p) corresponding to each lambda value.
# Process:
#   - Compute the correlation matrix using Corr.new2.
#   - Generate a sequence of lambda values from the largest eigenvalue of the correlation matrix to 0.
#   - For each lambda, apply ProxAlg_FGM to compute the support matrix.
# ------------------------------------------------------------
MGLASSO2 <- function(data, p, Ndim, M = 5, L = 100, ridge = ridge) {
  n <- dim(data)[1]
  tau <- dim(data)[2]
  est.cov.pc <- Corr.new2(xvec = data, ridge = ridge, Ndim = Ndim, p = p, M = 5)
  lambda.cmcg.max <- (eigen(est.cov.pc)$values)[1]
  lambda.cmcg <- seq(lambda.cmcg.max, 0, length.out = L)
  G.fglasso <- list()
  for (j in 1:L) {
    G.fglasso[[j]] <- ProxAlg_FGM(S = est.cov.pc, p = p, M = M, gamma = lambda.cmcg[j],
                                  Eta = "Auto", n.iteration = 10)$Support  # p x p support matrix
    print(j)
  }
  return(G.fglasso)
}

