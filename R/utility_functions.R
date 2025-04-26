# empirical measurements calculated from alpha, beta and gamma.
measure_emp <- function(alpha, beta, r){
  med.true <-  (beta != 0) & (alpha != 0)
  sig2_a.emp <- var(alpha[med.true])
  sig2_b.emp <- var(beta[med.true])*sum(med.true)
  varl.emp <- r^2 + sig2_a.emp*sig2_b.emp + 1
  Rmed.emp <- sig2_a.emp*sig2_b.emp/varl.emp
  Q.emp <- Rmed.emp/(Rmed.emp + r^2/varl.emp)
  return(c('Rmed.emp' = Rmed.emp, 'Q.emp'= Q.emp))
}

# FDR control method
est_sig = function(p_values, method = "BH", cutoff = 0.01) {
  if (method == "storey") {
    if (!requireNamespace("qvalue", quietly = TRUE)) {
      stop("Package 'qvalue' is required for Storey's method.")
    }
    qobj <- qvalue::qvalue(p_values)
    adjusted_p_values <- qobj$qvalues
  } else {
    adjusted_p_values <- p.adjust(p_values, method = method)
  }
  
  significant <- adjusted_p_values < cutoff
  return(significant)
}

#est_sig_BH = function(p_values, cutoff = 0.01){
#  adjusted_p_values <- p.adjust(p_values, method = "BH")
#  significant <- adjusted_p_values < cutoff
#  return(significant)
#}

# removal of highly correlated mediators.
remove_highly_correlated <- function(data, threshold = 0.8) {
  corr_matrix <- cor(data)
  to_keep <- rep(TRUE, ncol(data))

  # Loop through each column in the correlation matrix
  for (i in 1:(ncol(corr_matrix) - 1)) {
    for (j in (i + 1):ncol(corr_matrix)) {
      # Check if correlation exceeds the threshold and the column has not been removed yet
      if (abs(corr_matrix[i, j]) > threshold && to_keep[j]) {
        to_keep[j] <- FALSE  # Mark this column for removal
      }
    }
  }
  return(to_keep)
}

phi <- function(t){
  dnorm(t)
}

population_sd <- function(column, control_ind, K){
  return(sqrt((1-K)^2*var(column[control_ind]) + K^2*var(column[-control_ind])))
}

standardization <- function(X, M, Y, K, covariates){
  P <- mean(Y)
  p <- ncol(M)
  w <- ifelse(Y == 1, 1, P * (1 - K) / K / (1 - P))

  # standardize X
  control_ind <- which(Y == 0)
  X.stand <- scale(X, center = (1-K)*mean(X[control_ind]) + K*mean(X[-control_ind]), scale = sqrt((1-K)^2*var(X[control_ind]) + K^2*var(X[-control_ind])))

  # standardize M
  residual.M <- apply(M, 2, function(Var) {residuals(lm(Var ~ X + covariates, weights = w))})
  residual.sd <- apply(residual.M, 2, population_sd, control_ind = control_ind, K = K)
  M.stand <- scale(M, center = FALSE, scale = residual.sd)
  return(list('X' = X.stand, 'M' = M.stand))
}
