###############################################################
#                       MAPO
# This script implements the MAPO algorithm for functional data.
# It generates basis functions using k_T, computes the L2 inner product
# via Simpson's rule, and handles both balanced and unbalanced cases.
# In both cases, the time set 'ttt' is obtained from the output of
# 'generate x-p30-n300.R'; in the balanced case, all rows of 'ttt' are identical,
# while in the unbalanced case, rows are randomly generated.
###############################################################

########################################
# Function: nzmat
# Purpose: Remove zero rows from a matrix.
# Input:
#   - mat: A matrix.
# Output:
#   - A matrix with all-zero rows removed.
########################################
nzmat <- function(mat) {
  zero_rows <- apply(mat, 1, function(row) all(row == 0))
  mat_filtered <- mat[!zero_rows, ]
  return(mat_filtered)
}

###############################################################
# Function: simpson
# Purpose: Generate Simpson's weights.
# Input:
#   - ns: Number of sample points (must be odd).
# Output:
#   - A vector of Simpson's weights.
###############################################################
simpson <- function(ns) {
  wei <- rep(0, ns)
  wei[1] <- 1
  wei[ns] <- 1
  for (i in 1:((ns - 1) / 2)) {
    wei[2 * i] <- 4
  }
  for (i in 1:((ns - 1) / 2 - 1)) {
    wei[2 * i + 1] <- 2
  }
  h <- 1 / (ns - 1)
  wei <- wei * (h / 3)
  return(wei)
}

###############################################################
# Function: onorm
# Purpose: Compute the operator norm of a matrix.
# Input:
#   - a: A square matrix.
# Output:
#   - The largest eigenvalue of (a + t(a)) / 2.
###############################################################
onorm <- function(a) {
  return(eigen(round((a + t(a)) / 2, 8))$values[1])
}

###############################################################
# Function: tr
# Purpose: Compute the trace of a matrix.
# Input:
#   - a: A square matrix.
# Output:
#   - The sum of the diagonal elements of a.
###############################################################
tr <- function(a) {
  return(sum(diag(a)))
}

###############################################################
# Function: matpower
# Purpose: Compute the matrix power (Moore-Penrose type), keeping only the available eigenvalues.
# Input:
#   - a: A square matrix.
#   - alpha: The power to which the matrix is raised.
# Output:
#   - The matrix a raised to the power alpha.
###############################################################
matpower <- function(a, alpha) {
  a <- (a + t(a)) / 2
  tmp <- eigen(a)
  if (length(tmp$values) == 1) {
    m <- tmp$vectors %*% (abs(tmp$values)^alpha) %*% t(tmp$vectors)
  } else {
    m <- tmp$vectors %*% diag(abs(tmp$values)^alpha) %*% t(tmp$vectors)
  }
  return(m)
}

###############################################################
# Function: mppower
# Purpose: Compute the Moore-Penrose pseudoinverse type power of a matrix.
#          Eigenvalues below the 'ignore' threshold are discarded.
# Input:
#   - matrix: A square matrix.
#   - power: The exponent to which the eigenvalues are raised.
#   - ignore: Threshold for ignoring small eigenvalues.
# Output:
#   - The matrix raised to the specified power with small eigenvalues ignored.
###############################################################
mppower <- function(matrix, power, ignore) {
  eig <- eigen(matrix)
  eval <- eig$values
  evec <- eig$vectors
  m <- length(eval[abs(eval) > ignore])
  if (m == 1) {
    tmp <- (eval[1:m]^power) * evec[, 1:m] %*% t(evec[, 1:m])
  } else {
    tmp <- evec[, 1:m] %*% diag(eval[1:m]^power) %*% t(evec[, 1:m])
  }
  return(tmp)
}

###############################################################
# Function: gramt
# Purpose: Compute the Gram matrix for the first-level Hilbert space
#          based on a given kernel and time sets.
# Inputs:
#   - tte: Time points for evaluation.
#   - tto: Time points for observations.
#   - kern: Kernel type ('gauss' or 'brown').
# Output:
#   - The Gram matrix (transposed) computed according to the specified kernel.
###############################################################
gramt <- function(tte, tto, kern) {
  ltte <- length(tte)
  ltto <- length(tto)
  if (kern == "gauss") {
    a1 <- matrix(tte^2, ltte, ltto)
    a2 <- tte %*% t(tto)
    a3 <- t(matrix(tto^2, ltto, ltte))
    a <- a1 - 2 * a2 + a3
    b1 <- matrix(tto^2, ltto, ltto)
    b2 <- tto %*% t(tto)
    b3 <- t(matrix(tto^2, ltto, ltto))
    b <- b1 - 2 * b2 + b3
    sigma <- sum(sqrt(b)) / (2 * choose(ltto, 2))
    gamma <- 1 / (2 * sigma^2)
    ktmat <- exp(-gamma * a)
  }
  if (kern == "brown") {
    arr <- array(0, c(ltte, ltto, 2))
    arr[, , 1] <- matrix(tte, ltte, ltto)
    arr[, , 2] <- t(matrix(tto, ltto, ltte))
    ktmat <- apply(arr, c(1, 2), min)
  }
  return(t(ktmat))
}

###############################################################
# Function: gcvt.new
# Purpose: Compute the Generalized Cross Validation (GCV) score for epsilon_T.
# Inputs:
#   - N: Vector of sample sizes for each group.
#   - xxx: 3D array of observed random functions (n x nt x p).
#   - ttt: Time points corresponding to the observations.
#   - et: Regularization parameter.
#   - kern: Kernel type ('brown' or 'gauss').
# Output:
#   - The GCV score.
###############################################################
gcvt.new <- function(N, xxx, ttt, et, kern) {
  n <- dim(xxx)[1]
  nt <- dim(xxx)[2]
  p <- dim(xxx)[3]
  nuset <- numeric()
  deset <- numeric()
  
  for (j in 1:p) {
    for (i in 1:((n * N[j]) / max(N))) {
      if (dim(xxx)[2] == length(ttt)) {
        kt <- gramt(ttt, ttt, kern)
      } else {
        kt <- gramt(ttt[i, ], ttt[i, ], kern)
      }
      scale <- onorm(kt)
      ktinv <- matpower(kt + et * scale * diag(nt), -1)
      nuset <- c(nuset, sum((xxx[i, , j] - kt %*% ktinv %*% xxx[i, , j])^2))
      deset <- c(deset, (1 - tr(kt %*% ktinv) / nt)^2)
    }
  }
  
  out <- sum(nuset / deset)
  return(out)
}

###############################################################
# Function: evalx
# Purpose: Estimate function values at evaluation points using kernel regression.
# Inputs:
#   - f: Vector of function values at observation points.
#   - tte: Time points for evaluation.
#   - tto: Time points for observations.
#   - ridge: Regularization parameter.
#   - kern: Kernel type ('brown' or 'gauss').
# Output:
#   - A vector of estimated function values at the evaluation points.
###############################################################
evalx <- function(f, tte, tto, ridge, kern) {
  kt <- gramt(tto, tto, kern)
  scale <- eigen(kt)$values[1]
  ktinv <- matpower(kt + scale * ridge * diag(nrow(kt)), -1)
  kt1 <- t(gramt(tte, tto, kern))
  out <- kt1 %*% ktinv %*% f
  return(c(out))
}

###############################################################
# Function: evalxmat
# Purpose: Estimate a sample of functions.
# Inputs:
#   - ff: Matrix of function values (rows correspond to samples).
#   - tte: Time points for evaluation.
#   - ttt: Matrix of time points for observations (each row for a sample).
#   - ridge: Regularization parameter.
#   - kern: Kernel type ('brown' or 'gauss').
# Output:
#   - A matrix where each row contains the estimated function values for a sample.
###############################################################
evalxmat <- function(ff, tte, ttt, ridge, kern) {
  n <- dim(ff)[1]
  ffcoo <- numeric()
  for (i in 1:n) {
    ffcoo <- rbind(ffcoo, evalx(ff[i, ], tte, ttt[i, ], ridge, kern))
  }
  return(ffcoo)
}

###############################################################
# Function: gramx.new
# Purpose: Compute the Gram matrix for the second RKHS for a sample of functions.
# Inputs:
#   - Ni: Group size for the current variable.
#   - x: Matrix of estimated functions.
#   - kt: Diagonal matrix of Simpson's weights.
# Output:
#   - A Gram matrix computed based on the L2 inner product.
###############################################################
gramx.new <- function(Ni, x, kt) {
  n <- dim(x)[1] / Ni
  k2 <- x %*% kt %*% t(x)
  k1 <- t(matrix(diag(k2), n * Ni, n * Ni))
  k3 <- t(k1)
  k <- k1 - 2 * k2 + k3
  
  sigma <- sum(sqrt(k)) / (2 * choose(n * Ni, 2))
  gamma <- 1 / (2 * sigma^2)
  
  K <- matrix(data = 0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- sum(k[((1 + (i - 1) * Ni):(i * Ni)), ((1 + (j - 1) * Ni):(j * Ni))]) / (Ni^2)
    }
  }
  
  K <- exp(-gamma * K)
  return(K)
}

###############################################################
# Function: gramx.new2
# Purpose: Alternative computation of the Gram matrix for the second RKHS.
# Inputs:
#   - x: 3D array of estimated functions.
#   - kt: Diagonal matrix of Simpson's weights.
# Output:
#   - A Gram matrix computed from the rearranged data.
###############################################################
gramx.new2 <- function(x, kt) {
  n <- dim(x)[1]
  tau <- dim(x)[2]
  Ni <- dim(x)[3]
  x <- matrix(aperm(x, c(1, 3, 2)), n * Ni, tau)
  
  k2 <- x %*% kt %*% t(x)
  k1 <- t(matrix(diag(k2), dim(x)[1], dim(x)[1]))
  k3 <- t(k1)
  k <- k1 - 2 * k2 + k3
  
  sigma <- sum(sqrt(abs(k))) / (2 * choose(dim(x)[1], 2))
  gamma <- 1 / (2 * sigma^2)
  
  K <- matrix(data = 0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- sum(k[((1 + (i - 1) * Ni):(i * Ni)), ((1 + (j - 1) * Ni):(j * Ni))]) / (Ni^2)
    }
  }
  
  K <- exp(-gamma * K)
  return(K)
}

###############################################################
# Function: qmat
# Purpose: Compute the centering matrix Q = I - 1/n * 11'.
# Input:
#   - n: Number of observations.
# Output:
#   - The centering matrix Q.
###############################################################
qmat <- function(n) {
  return(diag(n) - (1/n) * (rep(1, n) %*% t(rep(1, n))))
}

###############################################################
# Function: cgram
# Purpose: Compute the centered Gram matrix G = Q*K*Q.
# Input:
#   - K: Original Gram matrix.
# Output:
#   - The centered Gram matrix.
###############################################################
cgram <- function(K) {
  n <- dim(K)[1]
  Q <- qmat(n)
  return(Q %*% K %*% Q)
}

###############################################################
# Function: aili
# Purpose: Compute matrices A_i and L_i for a given centered Gram matrix.
# Inputs:
#   - Gi: Centered Gram matrix for group i.
#   - ridge: Regularization parameter.
# Output:
#   - A list containing:
#       - Ai: Matrix A_i.
#       - Li: Matrix L_i.
###############################################################
aili <- function(Gi, ridge) {
  n <- dim(Gi)[1]
  scale <- eigen(Gi)$values[1]
  Q <- qmat(n)
  mat <- Gi + ridge * scale * Q
  Ai <- (n^(-1/2)) * mppower(Gi, 1/2, 1e-7) %*% mppower(mat, -1/2, 1e-7)
  Li <- Q - Ai %*% Ai
  return(list(Ai = Ai, Li = Li))
}

###############################################################
# Function: fapo.new
# Purpose: Compute the operator norms between groups using Simpson's rule
#          to approximate the L2 inner product (FAPO method).
# Inputs:
#   - N: Vector of sample sizes for each group.
#   - xvec: 3D array of estimated functions with dimensions (n*N) x tau x p.
#   - ridge: Regularization parameter.
# Output:
#   - A matrix of operator norms between groups.
###############################################################
fapo.new <- function(N, xvec, ridge) {
  n <- dim(xvec)[1] / max(N)
  p <- dim(xvec)[3]
  eps <- 1e-7
  ns <- dim(xvec)[2]
  kt <- diag(simpson(ns))
  H <- numeric()
  Linv <- numeric()
  
  for (i in 1:p) {
    Gi <- cgram(gramx.new(N[i], nzmat(xvec[, , i]), kt))
    store <- aili(Gi, ridge)
    H <- rbind(H, store$Ai)
    Linv <- rbind(Linv, mppower(store$Li, -1, eps))
  }
  
  M <- 0
  for (i in 1:p) {
    Hi <- H[((i - 1) * n + 1):(i * n), ]
    Linvi <- Linv[((i - 1) * n + 1):(i * n), ]
    M <- M + Hi %*% Linvi %*% Hi
  }
  
  Q <- qmat(n)
  Minv <- mppower(M + Q, -1, eps)
  LinvH <- numeric()
  for (i in 1:p) {
    LinvH <- rbind(LinvH, Linv[((i - 1) * n + 1):(i * n), ] %*% H[((i - 1) * n + 1):(i * n), ])
  }
  big <- LinvH %*% Minv %*% t(LinvH)
  norm1 <- matrix(0, p, p)
  for (i in 2:p) {
    for (j in 1:(i - 1)) {
      norm1[i, j] <- norm1[j, i] <- onorm(big[((i - 1) * n + 1):(i * n), ((j - 1) * n + 1):(j * n)])
    }
  }
  return(norm1)
}

###############################################################
# Function: fapo.new2
# Purpose: Alternative computation of the FAPO method for grouped functional data.
# Inputs:
#   - xvec: 3D array of estimated functions with dimensions n x tau x (N1+...+Np).
#   - ridge: Regularization parameter.
#   - Ndim: Vector of sample sizes for each group.
#   - p: Number of groups.
# Output:
#   - A matrix of operator norms between groups.
###############################################################
fapo.new2 <- function(xvec, ridge, Ndim, p) {
  n <- dim(xvec)[1]
  ns <- dim(xvec)[2]
  eps <- 1e-7
  kt <- diag(simpson(ns))
  H <- numeric()
  Linv <- numeric()
  
  for (i in 1:p) {
    index <- (sum(Ndim[1:(i - 1)]) + 1):(sum(Ndim[1:i]))
    Gi <- cgram(gramx.new2(xvec[, , index], kt))
    store <- aili(Gi, ridge)
    H <- rbind(H, store$Ai)
    Linv <- rbind(Linv, mppower(store$Li, -1, eps))
  }
  
  M <- 0
  for (i in 1:p) {
    Hi <- H[((i - 1) * n + 1):(i * n), ]
    Linvi <- Linv[((i - 1) * n + 1):(i * n), ]
    M <- M + Hi %*% Linvi %*% Hi
  }
  
  Q <- qmat(n)
  Minv <- mppower(M + Q, -1, eps)
  LinvH <- numeric()
  for (i in 1:p) {
    LinvH <- rbind(LinvH, Linv[((i - 1) * n + 1):(i * n), ] %*% H[((i - 1) * n + 1):(i * n), ])
  }
  big <- LinvH %*% Minv %*% t(LinvH)
  norm1 <- matrix(0, p, p)
  for (i in 2:p) {
    for (j in 1:(i - 1)) {
      norm1[i, j] <- norm1[j, i] <- as.numeric(norm(big[((i - 1) * n + 1):(i * n), ((j - 1) * n + 1):(j * n)], "F"))
    }
  }
  return(norm1)
}

###############################################################
# Function: gcvx.new
# Purpose: Compute the Generalized Cross Validation (GCV) score for epsilon_X.
# Inputs:
#   - N: Vector of sample sizes for each group.
#   - xxxeva: 3D array of estimated functions with dimensions (n*N) x tau x p.
#   - ex: Regularization parameter.
# Output:
#   - The GCV score.
###############################################################
gcvx.new <- function(N, xxxeva, ex) {
  ns <- dim(xxxeva)[2]
  n <- dim(xxxeva)[1] / max(N)
  p <- dim(xxxeva)[3]
  kt <- diag(simpson(ns))
  nuset <- numeric()
  deset <- numeric()
  
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      Gi <- gramx.new(N[i], nzmat(xxxeva[, , i]), kt)
      Gj <- gramx.new(N[j], nzmat(xxxeva[, , j]), kt)
      scale <- onorm(Gi)
      Gin2 <- matpower(Gi + ex * scale * diag(n), -2)
      nuset <- c(nuset, sum((Gi - Gj %*% Gi %*% Gin2 %*% Gi)^2))
      deset <- c(deset, (1 - tr(Gi %*% Gin2 %*% Gi) / n)^2)
    }
  }
  
  return(sum(nuset / deset))
}

###############################################################
# Function: gcvx.new2
# Purpose: Alternative GCV computation for epsilon_X handling different tensor shapes.
# Inputs:
#   - xxxeva: 3D array of estimated functions with dimensions n x tau x (N1+...+Np).
#   - ex: Regularization parameter.
#   - Ndim: Vector of sample sizes for each group.
# Output:
#   - The GCV score.
###############################################################
gcvx.new2 <- function(xxxeva, ex, Ndim) {
  ns <- dim(xxxeva)[2]
  n <- dim(xxxeva)[1]
  p <- length(Ndim)
  kt <- diag(simpson(ns))
  nuset <- numeric()
  deset <- numeric()
  
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      index1 <- (sum(Ndim[1:(i - 1)]) + 1):(sum(Ndim[1:i]))
      Gi <- cgram(gramx.new2(xxxeva[, , index1], kt))
      index2 <- (sum(Ndim[1:(j - 1)]) + 1):(sum(Ndim[1:j]))
      Gj <- cgram(gramx.new2(xxxeva[, , index2], kt))
      scale <- onorm(Gi)
      Gin2 <- matpower(Gi + ex * scale * diag(n), -2)
      nuset <- c(nuset, sum((Gi - Gj %*% Gi %*% Gin2 %*% Gi)^2))
      deset <- c(deset, (1 - tr(Gi %*% Gin2 %*% Gi) / n)^2)
    }
  }
  
  return(sum(nuset / deset))
}
