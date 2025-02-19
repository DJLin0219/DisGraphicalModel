###############################################################
#                     FARO
# This script implements the Functional Adaptive Regression for Ordination (FARO)
# and its parallel version (FAPO_parallel). It estimates graph structures from
# functional data via an ADMM-based group lasso method. The script computes
# kernel representations, selects the optimal regularization parameter, and
# performs neighborhood selection along a regularization path.
###############################################################

# Load necessary libraries for advanced matrix operations and parallel computing
library(Matrix)
library(parallel)

# ------------------------------------------------------------
# Function: lambda.sup
# Purpose: Compute the largest singular value of the product A.X' and A.Y.
#          This value is used to determine the maximum lambda for regularization.
# Inputs:
#   A.X : Matrix of predictor variables.
#   A.Y : Matrix of response variables.
# Output:
#   Numeric value representing the largest singular value.
# ------------------------------------------------------------
lambda.sup <- function(A.X, A.Y) {
  n <- nrow(A.X)          # Number of observations
  M <- ncol(A.Y)          # Number of response variables
  p <- ncol(A.X) / M      # Number of predictor groups
  candidates <- numeric(p)  # Vector to store candidate values
  
  for (j in 1:p) {
    A.X.j <- A.X[, ((j - 1) * M + 1):(j * M)]
    candidates[j] <- norm(t(A.X.j) %*% A.Y, "F") / n
  }
  
  return(max(candidates))
}

# ------------------------------------------------------------
# Function: RedKer2
# Purpose: Perform kernel reduction via eigenvalue decomposition.
#          The function reduces the dimensionality of a kernel matrix.
# Inputs:
#   K : Kernel matrix to be reduced.
#   C : Integer; number of top eigenvalues/eigenvectors to retain.
# Output:
#   A list containing:
#     SCORE: The reduced representation of the kernel matrix.
#     d    : The largest eigenvalue.
# ------------------------------------------------------------
RedKer2 <- function(K, C) {
  eig <- eigen(K)
  U <- eig$vectors   # Eigenvectors of K
  D <- eig$values    # Eigenvalues of K
  
  # Sort eigenvalues in descending order
  sorted_indices <- order(D, decreasing = TRUE)
  D_sorted <- D[sorted_indices]
  U_sorted <- U[, sorted_indices]
  
  # Select top C eigencomponents and compute scores
  COEF <- U_sorted[, 1:C]
  SCORE <- COEF %*% diag(sqrt(D_sorted[1:C]))
  
  # Extract the largest eigenvalue
  d <- D_sorted[1]
  
  list(SCORE = SCORE, d = d)
}

# ------------------------------------------------------------
# Main FARO Function
# Purpose: Perform Functional Adaptive Regression for Ordination (FARO) to estimate
#          graph structures from functional data.
# Inputs:
#   p          : Integer; number of variables (will be updated from h).
#   M          : Integer; number of basis functions or components.
#   n          : Integer; number of observations or samples.
#   tau        : Integer; number of time points or discretization points.
#   thres.ctrl : Numeric; threshold control parameter for neighbor recognition.
#   tol.abs    : Numeric; absolute tolerance for ADMM convergence.
#   tol.rel    : Numeric; relative tolerance for ADMM convergence.
#   L          : Integer; number of lambda values in the regularization path.
#   h          : 3D array of observed data with dimensions (n, tau, p).
# Output:
#   A list of adjacency matrices (each p x p) representing the estimated graph structure
#   for each lambda value.
# ------------------------------------------------------------
FARO <- function(p, M, n, tau, thres.ctrl, tol.abs, tol.rel, L, h) {
  # Update p and tau based on the dimensions of h
  p <- dim(h)[3]
  tau <- dim(h)[2]
  
  # Initialize response and predictor matrices
  A.Y <- matrix(NA, n, M)
  A.X <- matrix(NA, n, (p - 1) * M)
  
  # Create a group structure matrix and compute its norms
  d.array <- matrix(1, nrow = p, ncol = (p - 1) * M)
  d <- vector("list", p)
  norm.adj <- numeric(p)
  for (k in 1:p) {
    d[[k]] <- d.array[k, ]
    norm.adj[k] <- norm(d[[k]], "2")
  }
  
  # List to store the estimated adjacency matrices
  G.list.gY <- vector("list", L)
  
  # Default warm start matrices for the ADMM algorithm
  P.def <- matrix(0, (p - 1) * M, M)
  Q.def <- matrix(0.1, (p - 1) * M, M)
  U.def <- matrix(0.01, (p - 1) * M, M)
  
  # Initialize matrix for storing neighborhood selection vectors
  V.array <- matrix(0, nrow = 0, ncol = p)
  KX <- matrix(0, n, p * M)      # Matrix to store kernel scores
  kt <- diag(simpson(tau))       # Diagonal matrix of Simpson's weights
  
  # Compute kernel matrices and reduce dimensionality for each variable
  for (j in 1:p) {
    Gj <- cgram(gramx(h[, , j], kt))
    KX[, ((j - 1) * M + 1):(j * M)] <- RedKer2(Gj, C = M)$SCORE
  }
  
  # Compute the maximum lambda for the regularization path for each variable
  l.max <- numeric(p)
  for (j in 1:p) {
    all_cols <- 1:(M * p)
    exclude_cols <- ((j - 1) * M + 1):(j * M)
    keep_cols <- setdiff(all_cols, exclude_cols)
    A.Y <- KX[, exclude_cols]
    A.X <- KX[, keep_cols]
    l.max[j] <- lambda.sup(A.X, A.Y)
  }
  lambda.max <- max(l.max)
  lambdas.gY <- seq(lambda.max, 0, length.out = L)
  
  # Loop over each variable to perform neighborhood selection along the lambda path
  for (j in 1:p) {
    all_cols <- 1:(M * p)
    exclude_cols <- ((j - 1) * M + 1):(j * M)
    keep_cols <- setdiff(all_cols, exclude_cols)
    A.Y <- KX[, exclude_cols]
    A.X <- KX[, keep_cols]
    
    P <- P.def
    Q <- Q.def
    U <- U.def
    V.j <- matrix(NA, L, p)
    
    for (l in 1:L) {
      lambda <- lambdas.gY[l]
      grp.lasso.result <- ADMM.grplasso(
        A.X = A.X, A.Y = A.Y, d = d[[j]],
        lambda = lambda, rho.init = 1,
        P.in = P, Q.in = Q, U.in = U,
        tol.abs = tol.abs, tol.rel = tol.rel,
        maxiter = 50
      )
      P <- grp.lasso.result$P
      Q <- grp.lasso.result$Q
      U <- grp.lasso.result$U
      
      # Process the estimated coefficient matrix P into a neighborhood selection vector
      P.frob <- numeric(p - 1)
      for (k in 1:(p - 1)) {
        P.frob[k] <- norm(P[((k - 1) * M + 1):(k * M), ], "F")
      }
      threshold <- lambda * thres.ctrl
      
      # Determine active neighbors based on the threshold
      V.jl <- numeric(p)
      for (juliet in 1:(p - 1)) {
        if (juliet < j) {
          if (P.frob[juliet] > threshold)
            V.jl[juliet] <- 1
        } else {
          if (P.frob[juliet] > threshold)
            V.jl[juliet + 1] <- 1
        }
      }
      V.j[l, ] <- V.jl
    }
    V.array <- rbind(V.array, V.j)
  }
  
  # Assemble the list of adjacency matrices from the neighborhood selection vectors
  G.list.gY <- vector("list", L)
  for (l in 1:L) {
    G.list.gY[[l]] <- matrix(NA, p, p)
    for (j in 1:p) {
      G.list.gY[[l]][j, ] <- V.array[((j - 1) * L + l), ]
    }
  }
  return(G.list.gY)
}

# ------------------------------------------------------------
# Function: FAPO_parallel
# Purpose: Parallel version of the FARO function. This function computes the neighborhood
#          selection vectors for a given lambda index across variables.
# Inputs:
#   p          : Integer; number of variables (updated from h).
#   M          : Integer; number of basis functions or components.
#   n          : Integer; number of observations.
#   tau        : Integer; number of time points.
#   thres.ctrl : Numeric; threshold control parameter.
#   tol.abs    : Numeric; absolute tolerance for ADMM convergence.
#   tol.rel    : Numeric; relative tolerance for ADMM convergence.
#   l          : Integer; index of the current lambda value in the regularization path.
#   L          : Integer; total number of lambda values (default is 100).
#   h          : 3D array of observed data with dimensions (n, tau, p).
# Output:
#   Matrix representing the neighborhood selection vector for each variable at the specified lambda.
# ------------------------------------------------------------
FAPO_parallel <- function(p, M, n, tau, thres.ctrl, tol.abs, tol.rel, l, L = 100, h) {
  # Update p and tau based on h's dimensions
  p <- dim(h)[3]
  tau <- dim(h)[2]
  
  # Initialize matrices and variables
  A.Y <- matrix(NA, n, M)
  A.X <- matrix(NA, n, (p - 1) * M)
  d.array <- matrix(1, nrow = p, ncol = (p - 1) * M)
  d <- vector("list", p)
  norm.adj <- numeric(p)
  for (k in 1:p) {
    d[[k]] <- d.array[k, ]
    norm.adj[k] <- norm(d[[k]], "2")
  }
  
  # Default warm start matrices for ADMM
  P.def <- matrix(0, (p - 1) * M, M)
  Q.def <- matrix(0.1, (p - 1) * M, M)
  U.def <- matrix(0.01, (p - 1) * M, M)
  
  # Initialize matrix for neighborhood selection vectors and compute kernel scores
  V.array <- matrix(0, nrow = 0, ncol = p)
  KX <- matrix(0, n, p * M)
  kt <- diag(simpson(tau))
  for (j in 1:p) {
    Gj <- cgram(gramx(h[, , j], kt))
    KX[, ((j - 1) * M + 1):(j * M)] <- RedKer2(Gj, C = M)$SCORE
  }
  
  # Compute maximum lambda for each variable and generate the lambda sequence
  l.max <- numeric(p)
  for (j in 1:p) {
    all_cols <- 1:(M * p)
    exclude_cols <- ((j - 1) * M + 1):(j * M)
    keep_cols <- setdiff(all_cols, exclude_cols)
    A.Y <- KX[, exclude_cols]
    A.X <- KX[, keep_cols]
    l.max[j] <- lambda.sup(A.X, A.Y)
  }
  lambda.max <- max(l.max)
  lambdas.gY <- seq(lambda.max, 0, length.out = L)
  
  # For each variable, compute the neighborhood selection vector for the given lambda index
  for (j in 1:p) {
    all_cols <- 1:(M * p)
    exclude_cols <- ((j - 1) * M + 1):(j * M)
    keep_cols <- setdiff(all_cols, exclude_cols)
    A.Y <- KX[, exclude_cols]
    A.X <- KX[, keep_cols]
    
    P <- P.def
    Q <- Q.def
    U <- U.def
    V.j <- matrix(NA, 1, p)
    
    lambda <- lambdas.gY[l]
    grp.lasso.result <- ADMM.grplasso(
      A.X = A.X, A.Y = A.Y, d = d[[j]],
      lambda = lambda, rho.init = 1,
      P.in = P, Q.in = Q, U.in = U,
      tol.abs = tol.abs, tol.rel = tol.rel,
      maxiter = 50
    )
    P <- grp.lasso.result$P
    Q <- grp.lasso.result$Q
    U <- grp.lasso.result$U
    
    # Convert estimated P into a neighborhood selection vector
    P.frob <- rep(0, p - 1)
    for (k in 1:(p - 1)) {
      P.frob[k] <- norm(P[(k - 1) * M + (1:M), ], "F")
    }
    threshold <- lambda * thres.ctrl
    V.jl <- rep(0, p)
    for (juliet in 1:(p - 1)) {
      if (juliet < j) {
        if (P.frob[juliet] > threshold)
          V.jl[juliet] <- 1
      } else {
        if (P.frob[juliet] > threshold)
          V.jl[juliet + 1] <- 1
      }
    }
    V.array <- rbind(V.array, V.jl)
  }
  
  return(V.array)
}
