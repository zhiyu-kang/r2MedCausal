#' @title Generate Synthetic Data for High-Dimensional Mediation Analysis
#' @description Simulates exposure, mediators, and binary outcomes with specified variance components. The function generates mediator effects and performs case-control sampling.
#'
#' @param p Integer. Number of mediators.
#' @param N Integer. Total sample size.
#' @param K Numeric. Proportion of cases in the population.
#' @param r Numeric. Direct effect of exposure on the outcome.
#' @param res.var Numeric. Residual variance in the outcome.
#' @param pi_11 Numeric. Proportion of mediators with a nonzero effect for both alpha and beta.
#' @param pi_alpha Numeric. Proportion of mediators with a nonzero alpha effect.
#' @param pi_beta Numeric. Proportion of mediators with a nonzero beta effect.
#' @param sig2_a Numeric. Variance of alpha effects.
#' @param sig2_11 Numeric. Variance of beta effects for shared mediators.
#' @param sig2_01 Numeric. Variance of beta effects for independent mediators.
#' @param alpha Numeric vector. Directly specified alpha values (optional).
#' @param beta Numeric vector. Directly specified beta values (optional).
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{x}{Numeric vector. Exposure values.}
#'   \item{y}{Binary outcome vector (0/1).}
#'   \item{liab}{Numeric vector. Underlying liability scores.}
#'   \item{M}{Matrix. Mediator values.}
#'   \item{Q}{Q measurement.}
#'   \item{Rmed}{Estimated mediation effect (theoretical).}
#'   \item{varl}{Variance of liability.}
#'   \item{Q.emp}{Empirical Q measurement.}
#'   \item{Rmed.emp}{Empirical mediation effect.}
#'   \item{alpha}{Vector of generated alpha values.}
#'   \item{beta}{Vector of generated beta values.}
#' }
#'
#' @examples
#' data <- data_generate(p = 1000, N = 100, K = 0.05, r = 0.2, res.var = 0.5,
#'                       pi_11 = 0.09, pi_alpha = 0.3, pi_beta = 0.3,
#'                       sig2_a = 0.5, sig2_11 = 0.2, sig2_01 = 0.3, seed = 123)
#' str(data)
#'
#' @export

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
  S_chol <- generate_S_chol(p, alpha, seed)
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
  true_emp <- measure_emp(alpha, beta, r, res.var)
  Rmed.emp <- unname(true_emp['Rmed.emp'])
  Q.emp <- unname(true_emp['Q.emp'])
  varl.emp <- unname(true_emp['varl.emp'])

  # Compute theoretical mediation effects only if sig2_a, sig2_11, sig2_01 were provided
  if (!is.null(sig2_a) && !is.null(sig2_11) && !is.null(sig2_01)) {
    varl <- r^2 + sig2_a * sig2_11 + sig2_11 + sig2_01 + res.var
    Rmed <- sig2_a * sig2_11 / varl
    Q <- Rmed / (Rmed + r^2 / varl)
  } else {
    varl <- NA
    Rmed <- NA
    Q <- NA
  }

  return(list('x' = x, 'y' = Yb, 'liab' = Y, 'M' = M, 'Q' = Q, 'Rmed' = Rmed, 'varl' = varl,
              'Q.emp' = Q.emp, 'Rmed.emp' = Rmed.emp, 'varl.emp' = varl.emp, 'alpha' = alpha, 'beta' = beta))
}

#' @keywords internal
generate_S_chol = function(p, alpha){
  set.seed(seed)
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

#' @keywords internal
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
