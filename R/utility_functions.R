
measure_emp <- function(alpha, beta, r){
  # empirical variances calculated from alpha and beta
  med.true <-  (beta != 0) & (alpha != 0)
  sig2_a.emp <- var(alpha[med.true])
  sig2_b.emp <- var(beta[med.true])*sum(med.true)
  varl.emp <- r^2 + sig2_a.emp*sig2_b.emp + 1
  Rmed.emp <- sig2_a.emp*sig2_b.emp/varl.emp
  V.emp <- Rmed.emp/(Rmed.emp + r^2/varl.emp)
  return(c('Rmed.emp' = Rmed.emp, 'V.emp'= V.emp))
}
