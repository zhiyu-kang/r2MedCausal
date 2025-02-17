# need to add Roxygen2 documentation

library(dplyr)

generate_S_chol = function(p, alpha){
  #var-cov of M and x
  M_cov <- matrix(0, p, p)
  xi_sd <- rep(1, p) # sd of uncorrelated components xi
  X_sd <- 1 # sd of X
  k <- 2 # num of latent factors (k >= 1)
  rho <- 0.5 # Should be less than 1;
  a_sd <- 0.3 # alpha sd

  if (k > 1) { # k > 1 means there is confounding in xi
    A.U <- matrix(rnorm(p * (k-1), sd = a_sd), nrow = p) #coefficient of U, varU=1
    M_cov <- alpha %*% t(alpha) * X_sd^2 / rho + A.U %*% t(A.U) + diag(xi_sd)
  } else {
    M_cov <- alpha %*% t(alpha) * X_sd^2 / rho + diag(xi_sd)
  }

  S <- matrix(NA, p+1, p+1)
  S[1,1] <- X_sd^2
  S[1,-1] <- alpha * X_sd^2
  S[-1,1] <- alpha * X_sd^2
  S[-1,-1] <- M_cov

  # Cache this matrix
  S_chol <- chol(S)
  return(S_chol)
}

generate_nonzero_ind <- function(pi_11, pi_a, pi_b, p, seed){
  set.seed(seed)
  n_a <- round(p * pi_a)
  n_b <- round(p * pi_b)
  n_11 <- round(p * pi_11)

  ind_11 = sort(sample(p, n_11, replace = F))
  ind_10 = sort(sample(setdiff(seq_len(p), ind_11), n_a - n_11, replace = F))
  ind_a = sort(c(ind_11, ind_10))
  ind_01 = sort(sample(setdiff(seq_len(p), ind_a), n_b - n_11, replace = F))
  ind_b = sort(c(ind_11, ind_01))

  return(list(ind_a = ind_a, ind_b = ind_b, ind_11 = ind_11, ind_01 = ind_01, ind_10 = ind_10))
}

# main function for data generation
data_generate <- function(p, N, K, r, res.var,
                          pi_11 = NULL, pi_alpha = NULL, pi_beta = NULL,
                          sig2_a = NULL, sig2_11 = NULL, sig2_01 = NULL,
                          alpha = NULL, beta = NULL, seed=123) {

  set.seed(seed)

  # Check if alpha and beta are provided
  if (is.null(alpha) || is.null(beta)) {
    if (is.null(sig2_a) || is.null(sig2_11) || is.null(sig2_01) || is.null(pi_11) || is.null(pi_alpha) || is.null(pi_beta)) {
      stop("If alpha and beta are not provided, sig2_a, sig2_11, sig2_01 pi_11, pi_alpha, and pi_beta must be specified.")
    }

    # Generate nonzero indices
    nonzero_ind <- generate_nonzero_ind(pi_11, pi_alpha, pi_beta, p, seed)

    # Generate alpha
    alpha <- rnorm(p, 0, sqrt(sig2_a))
    alpha[-nonzero_ind$ind_a] <- 0

    # Generate beta
    beta <- rep(0, p)
    beta[nonzero_ind$ind_11] <- rnorm(length(nonzero_ind$ind_11), 0, sqrt(sig2_11 / length(nonzero_ind$ind_11)))
    beta[nonzero_ind$ind_01] <- rnorm(length(nonzero_ind$ind_01), 0, sqrt(sig2_01 / length(nonzero_ind$ind_01)))
  }

  # Generate M and X using correlation structure
  S_chol <- generate_S_chol(p, alpha)
  XM <- matrix(rnorm(N/K/2 * (p+1)), nrow = N/K/2)
  XM <- XM %*% S_chol
  x <- XM[,1]
  M <- XM[,-1]

  # Generate Y (liability) and Yb (binary outcome)
  Y <- r*x + as.vector(M %*% matrix(beta)) + rnorm(N/K/2, 0, sqrt(res.var))
  varl <- var(Y)
  Yb <- ifelse(Y > qnorm(1-K, 0, sqrt(varl)), 1, 0)

  # Case-control sampling with P=0.5
  if (sum(Yb == 1) > N/2) {
    caseind <- sample(which(Yb == 1), N/2)
  } else {
    caseind <- which(Yb == 1)
  }

  controlind <- which(Yb == 0)
  if (length(controlind) > (N - length(caseind))) {
    controlind <- sample(controlind, N - length(caseind))
  }

  selectind <- sort(c(caseind, controlind))
  x <- x[selectind]
  M <- M[selectind, ]
  Yb <- Yb[selectind]
  Y <- Y[selectind]

  # Compute empirical mediation effects
  true_emp <- measure_emp(alpha, beta, r)
  Rmed.emp <- unname(true_emp['Rmed.emp'])
  V.emp <- unname(true_emp['V.emp'])

  # Compute theoretical mediation effects only if sig2_a, sig2_11, sig2_01 were provided
  if (!is.null(sig2_a) && !is.null(sig2_11) && !is.null(sig2_01)) {
    varl <- r^2 + sig2_a * sig2_11 + sig2_11 + sig2_01 + res.var
    Rmed <- sig2_a * sig2_11 / varl
    V <- Rmed / (Rmed + r^2 / varl)
  } else {
    varl <- NA
    Rmed <- NA
    V <- NA
  }

  return(list('x' = x, 'y' = Yb, 'liab' = Y, 'M' = M, 'V' = V, 'Rmed' = Rmed, 'varl' = varl,
              'V.emp' = V.emp, 'Rmed.emp' = Rmed.emp, 'alpha' = alpha, 'beta' = beta))
}
