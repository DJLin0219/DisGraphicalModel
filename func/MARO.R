###############################################################
#                      [Date: YYYY/MM/DD]   MARO
# This script implements the MARO algorithm for estimating
# graph structures from functional data via ADMM-based group lasso.
# It includes functions for computing singular values, kernel reduction,
# and two main routines (MARO and MAPO2) for graph estimation.
###############################################################

# Load necessary libraries
library(Matrix)    # For matrix operations
library(parallel)  # For parallel computing

# ------------------------------------------------------------
# Function: lambda.sup
# Purpose: Compute the largest singular value of the product of A.X' and A.Y.
# Inputs:
#   A.X : Matrix of predictor variables.
#   A.Y : Matrix of response variables.
# Output:
#   The largest singular value.
# ------------------------------------------------------------
lambda.sup <- function(A.X, A.Y) {
  n <- nrow(A.X)  # Number of observations
  M <- ncol(A.Y)  # Number of response variables
  p <- ncol(A.X) / M  # Number of predictor groups
  candidates <- numeric(p)  # Initialize vector to store candidate values
  
  for (j in 1:p) {
    A.X.j <- A.X[, ((j - 1) * M + 1):(j * M)]
    candidates[j] <- norm(t(A.X.j) %*% A.Y, "F") / n
  }
  
  return(max(candidates))
}

# ------------------------------------------------------------
# Function: RedKer2
# Purpose: Perform kernel reduction using eigenvalue decomposition.
# Inputs:
#   K : Kernel matrix to be reduced.
#   C : Number of top eigenvalues and corresponding eigenvectors to retain.
# Output:
#   A list containing:
#       SCORE: Reduced representation of the kernel matrix.
#       d    : The largest eigenvalue.
# ------------------------------------------------------------
RedKer2 <- function(K, C) {
  eig <- eigen(K)
  U <- eig$vectors   # Eigenvectors
  D <- eig$values    # Eigenvalues
  
  # Sort eigenvalues in descending order
  sorted_indices <- order(D, decreasing = TRUE)
  D_sorted <- D[sorted_indices]
  U_sorted <- U[, sorted_indices]
  
  # Select top C eigenvalues and corresponding eigenvectors
  COEF <- U_sorted[, 1:C]
  SCORE <- COEF %*% diag(sqrt(D_sorted[1:C]))
  
  # Extract the largest eigenvalue
  d <- D_sorted[1]
  
  list(SCORE = SCORE, d = d)
}

# ------------------------------------------------------------
# Function: MARO
# Purpose: Main MAPO function to estimate graph structure when the data
#          has dimensions (n*N) x tau x p.
# Description:
#   This function processes a 3-dimensional array of functional data with
#   dimensions (n, tau, p) and, for each variable, estimates its neighborhood
#   structure over a regularization path defined by L lambda values.
# Inputs:
#   p          : Integer. Number of variables (dimensionality; updated from h).
#   M          : Integer. Number of basis functions or components.
#   n          : Integer. Number of observations or samples.
#   tau        : Integer. Number of time points (discretization points).
#   thres.ctrl : Numeric. Threshold control parameter for neighbor recognition.
#   tol.abs    : Numeric. Absolute tolerance for convergence in the ADMM algorithm.
#   tol.rel    : Numeric. Relative tolerance for convergence in the ADMM algorithm.
#   L          : Integer. Number of lambda values in the regularization path.
#   h          : Array. 3-dimensional array of observed data with dimensions (n, tau, p).
# Output:
#   A list of adjacency matrices representing the estimated graph structure
#   for each lambda value.
# ------------------------------------------------------------
MARO <- function(p, M, n, tau, thres.ctrl, tol.abs, tol.rel, L, h) {
  p <- dim(h)[3]   # Update p based on the third dimension of h
  tau <- dim(h)[2] # Update tau based on the second dimension of h
  A.Y <- matrix(NA, n, M)           # Initialize response matrix
  A.X <- matrix(NA, n, (p - 1) * M)   # Initialize predictor matrix
  d.array <- matrix(1, nrow = p, ncol = (p - 1) * M)  # Group structure matrix
  d <- vector("list", p)            # List to store group structures
  norm.adj <- numeric(p)            # Vector to store group norms
  
  for (k in 1:p) {
    d[[k]] <- d.array[k, ]
    norm.adj[k] <- norm(d[[k]], "2")
  }
  
  G.list.gY <- vector("list", L)    # List to store adjacency matrices
  
  # Default warm start matrices for ADMM algorithm
  P.def <- matrix(0, (p - 1) * M, M)
  Q.def <- matrix(0.1, (p - 1) * M, M)
  U.def <- matrix(0.01, (p - 1) * M, M)
  
  # Initialize matrix to store neighborhood selection vectors
  V.array <- matrix(0, nrow = 0, ncol = p)
  KX <- matrix(0, n, p * M)         # Matrix to store kernel scores
  kt <- diag(simpson(tau))          # Diagonal matrix of Simpson's weights
  
  # Compute kernel matrices and reduce dimensionality for each variable
  for (j in 1:p) {
    Gj <- cgram(gramx(h[, , j], kt))
    KX[, ((j - 1) * M + 1):(j * M)] <- RedKer2(Gj, C = M)$SCORE
  }
  
  # Compute maximum lambda for regularization path
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
  
  # Main loop: For each variable, perform neighborhood selection over lambda path
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
        P.frob[k] <- norm(P[(k - 1) * M + (1:M), ], "F")
      }
      
      threshold <- lambda * thres.ctrl
      
      # Neighbor recognition: determine active neighbors based on threshold
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
      
      V.j[l, ] <- V.jl
    }
    V.array <- rbind(V.array, V.j)
  }
  
  # Construct the list of adjacency matrices from the neighborhood vectors
  G.list.gY <- list()
  for (l in 1:L) {
    G.list.gY[[l]] <- matrix(NA, p, p)
    for (j in 1:p) {
      G.list.gY[[l]][j, ] <- V.array[(j - 1) * L + l, ]
    }
  }
  return(G.list.gY)
}

# ------------------------------------------------------------
# Function: MAPO2
# Purpose: Main MAPO2 function to estimate graph structure when the data
#          has dimensions (n) x tau x (N1 + ... + Np).
# Description:
#   This function processes a 3-dimensional array of observed data and applies
#   ADMM-based group lasso to estimate the graph structure for each lambda
#   along a regularization path.
# Inputs:
#   p          : Integer. Number of variables (dimensionality; updated from h).
#   M          : Integer. Number of basis functions or components.
#   n          : Integer. Number of observations or samples.
#   tau        : Integer. Number of time points (discretization points).
#   thres.ctrl : Numeric. Threshold control parameter for neighbor recognition.
#   tol.abs    : Numeric. Absolute tolerance for convergence in the ADMM algorithm.
#   tol.rel    : Numeric. Relative tolerance for convergence in the ADMM algorithm.
#   L          : Integer. Number of lambda values in the regularization path.
#   h          : Array. 3-dimensional array of observed data with dimensions (n, tau, p).
# Output:
#   A list of adjacency matrices representing the estimated graph structure
#   for each lambda value.
# ------------------------------------------------------------
MAPO2 <- function(p, M, n, tau, thres.ctrl, tol.abs, tol.rel, L, h) {
  p <- dim(h)[3]   # Update p based on the third dimension of h
  tau <- dim(h)[2] # Update tau based on the second dimension of h
  A.Y <- matrix(NA, n, M)           # Initialize response matrix
  A.X <- matrix(NA, n, (p - 1) * M)   # Initialize predictor matrix
  d.array <- matrix(1, nrow = p, ncol = (p - 1) * M)  # Group structure matrix
  d <- vector("list", p)            # List to store group structures
  norm.adj <- numeric(p)            # Vector to store group norms
  
  for (k in 1:p) {
    d[[k]] <- d.array[k, ]
    norm.adj[k] <- norm(d[[k]], "2")
  }
  
  G.list.gY <- vector("list", L)    # List to store adjacency matrices
  
  # Default warm start matrices for ADMM algorithm
  P.def <- matrix(0, (p - 1) * M, M)
  Q.def <- matrix(0.1, (p - 1) * M, M)
  U.def <- matrix(0.01, (p - 1) * M, M)
  
  # Initialize matrix to store neighborhood selection vectors
  V.array <- matrix(0, nrow = 0, ncol = p)
  KX <- matrix(0, n, p * M)         # Matrix to store kernel scores
  kt <- diag(simpson(tau))          # Diagonal matrix of Simpson's weights
  
  # Compute kernel matrices and reduce dimensionality for each variable
  for (j in 1:p) {
    Gj <- cgram(gramx(h[, , j], kt))
    KX[, ((j - 1) * M + 1):(j * M)] <- RedKer2(Gj, C = M)$SCORE
  }
  
  # Compute maximum lambda for regularization path
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
  
  # Main loop: For each variable, perform neighborhood selection over lambda path
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
        P.frob[k] <- norm(P[(k - 1) * M + (1:M), ], "F")
      }
      
      threshold <- lambda * thres.ctrl
      
      # Neighbor recognition: determine active neighbors based on threshold
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
      
      V.j[l, ] <- V.jl
    }
    V.array <- rbind(V.array, V.j)
  }
  
  # Construct the list of adjacency matrices from the neighborhood vectors
  G.list.gY <- list()
  for (l in 1:L) {
    G.list.gY[[l]] <- matrix(NA, p, p)
    for (j in 1:p) {
      G.list.gY[[l]][j, ] <- V.array[(j - 1) * L + l, ]
    }
  }
  return(G.list.gY)
}
