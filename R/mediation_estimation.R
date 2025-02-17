# change data generation process so that measurement can be compared directly.

######## add type 1,2 and 3 non-mediators #########
#library(ggplot2)
library(foreach)
library(doParallel)
#library(caret)
#library(MASS)
#library(broom)
#library(adaHuber)

##### estimation procedure #####

#################### adjust pcgc function and est_procedure function to incorporate covariates

# pcgc algorithm
pcgc = function(genotype.stand, phenotype, covariates=NA, P, K, P_cond=NA, adjustGRM = TRUE, Reg.model='ols'){
  phi=function(t){dnorm(t)}
  #calculate and adjust genetic correlation matrix
  if (adjustGRM == TRUE){
    X = as.matrix(covariates)
    genotype.stand = (diag(length(phenotype)) - X %*% solve(t(X) %*% X) %*% t(X)) %*% genotype.stand
  }
  correlation=genotype.stand %*%t (genotype.stand) / ncol(genotype.stand)
  w = ifelse(phenotype==1, 1, P*(1-K)/K/(1-P))

  data <- cbind(phenotype, covariates)
  logit = glm(data[, 1] ~ data[, -1], family = 'binomial')
  P_cond=predict(logit, type='response') # calculate "unbiased" P_cond first
  K_cond=(K*(1-P)*P_cond/(P*(1-K)))/(1+K*(1-P)*P_cond/(P*(1-K))-P_cond)
  t_cond=unname(qnorm(1-K_cond))
  # calculate sig2_t using law of total variance.
  VEt=K*(1-K)*(mean(t_cond[which(phenotype==1)]) - mean(t_cond[which(phenotype==0)]))^2
  EVt=K*var(t_cond[which(phenotype==1)]) + (1-K)*var(t_cond[which(phenotype==0)])
  sig2_t=VEt+EVt

  #phenotype correlation matrix
  product = matrix((phenotype-P_cond)/sqrt(P_cond*(1-P_cond)))%*%t(matrix((phenotype-P_cond)/sqrt(P_cond*(1-P_cond))))

  # calculate coefficient of regressing phenotypic correlation on genetypic correlation
  coef_vect=matrix(dnorm(t_cond)*(1-P_cond*(P-K)/(P*(1-K)))/sqrt(P_cond*(1-P_cond))/(K_cond+(1-K_cond)*K*(1-P)/P/(1-K)))
  coef_mat=coef_vect%*%t(coef_vect)
  predictor=coef_mat*correlation

  # regression using off-diagonal elements,  here Gij is different from the paper.
  Zij=c(product)[lower.tri(product)]; Gij=c(predictor)[lower.tri(product)]
  if (Reg.model=='ols'){
    model=lm(Zij~-1+Gij)
    sig2_g=coefficients(model)
  } else if (Reg.model=='adaHuber'){
    model = adaHuber.reg(as.matrix(Gij), Zij, method = "adaptive")
    sig2_g=model$coef[2]
  }
  return(list('sig2_g' = sig2_g, 'sig2_t' = sig2_t, 'P_cond' = P_cond, 'gamma' = gamma))
}

est_sig_BH = function(p_values, cutoff = 0.01){
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  significant <- adjusted_p_values < cutoff
  return(significant)
}

removeHighlyCorrelated <- function(data, threshold = 0.8) {
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

population_sd <- function(column, control_ind, K){
    return(sqrt((1-K)^2*var(column[control_ind]) + K^2*var(column[-control_ind])))
}

filter_type2_nonmediator <- function(X, M, Y, K, covariates_demo, W){
  P <- mean(Y)
  N <- length(Y)
  p <- ncol(M)
  w <- ifelse(Y == 1, 1, P * (1 - K) / K / (1 - P))
  covariates <- cbind(covariates_demo, W)

  residual.est <- apply(M, 2, function(Var) {residuals(lm(Var ~ X + covariates, weights = w))}) # new residual with mean 0 and sd 1 in population

  ind_cor <- removeHighlyCorrelated(residual.est, threshold = 0.8)

  getAlpha <- function(Var){
    model <- lm(Var ~ X + covariates_demo, weights = w)
    alpha.est <- summary(model)$coefficients[2,'Estimate']
    p_value <- summary(model)$coefficients[2,'Pr(>|t|)']
    return(list(alpha.est, p_value))
  }

  temp <- apply(M, 2, getAlpha)
  temp2 <- sapply(temp, function(x) sapply(x, function(y) y[[1]]))
  alpha.est <- temp2[1, ]
  p_value_alpha <- temp2[2, ]

  bias <- p*sum((w*X)^2)/(sum(w*(X)^2))^2
  alpha2.est <- sum(alpha.est^2) - bias # estimate of alpha^T alpha after bias correction

  #q_value_alpha <- qvalue(p_value_alpha)
  sig_alpha <- est_sig_BH(p_value_alpha)
  sig2_a.est <- alpha2.est/p/mean(sig_alpha)

  return(list('selected_var' = unname(sig_alpha & ind_cor), 'sig2_a.est' = sig2_a.est))
}

estimation <- function(X, M, Y, K, covariates_demo, W, filter_result) {
  P <- mean(Y)
  N <- length(Y)
  p <- ncol(M)
  w <- ifelse(Y == 1, 1, P * (1 - K) / K / (1 - P))
  covariates <- cbind(covariates_demo, W)

  residual.est <- apply(M, 2, function(Var) {residuals(lm(Var ~ X + covariates, weights = w))}) # new residual with mean 0 and sd 1 in population

  # filter out non-mediators
  sig_alpha <- filter_result$selected_var
  sig2_a.est <- filter_result$sig2_a.est
  residual.est.filter <- residual.est[, sig_alpha]

  # estimate sig2_11

  pcgc.est <- pcgc(residual.est.filter, Y, covariates = cbind(X, covariates), P = P, K = K, P_cond = NA, adjustGRM=TRUE, Reg.model='ols')
  sig2_11.est <- pcgc.est$sig2_g

  # fixed effect glm
  probit.fit <- glm(Y ~ X + covariates, family = binomial(link='probit'), weights = w)
  r.est <- probit.fit$coefficients[2]

  # estimate of Rmed and total variance
  varl.est <- r.est^2 + sig2_a.est*sig2_11.est + 1
  Rmed.est <- sig2_a.est*sig2_11.est/varl.est
  V.est <- Rmed.est/(Rmed.est + r.est^2/varl.est)

  return(c('Rmed.est' = unname(Rmed.est), 'V.est' = unname(V.est), 'varl.est' = unname(varl.est),
            'sig2_11.est' = unname(sig2_11.est), 'sig2_a.est' = unname(sig2_a.est), 'gamma.est' = unname(r.est),
            'sig_number' = sum(sig_alpha)))
}

normal_fit <- function(X, M, Y, K, covariates_demo) {
  start.time <- Sys.time()
  A <- scale(M) %*% t(scale(M))/ncol(M)
  eig_A <- eigen(A, symmetric = TRUE)
  U <- eig_A$vectors
  W <- sqrt(nrow(M)) * U[, 1:3]
  filter_result <- filter_type2_nonmediator(X, M, Y, K, covariates_demo, W)
  result <- estimation(X, M, Y, K, covariates_demo, W, filter_result)
  time.total <- as.numeric(Sys.time() - start.time, units = "secs")
  result <- c(result, 'time' = time.total)
  return(result)
}

cross_fit <- function(X, M, Y, K, covariates_demo, nfold = 5, seed = 123) {
  start.time <- Sys.time()
  set.seed(seed)
  ind <- sample(length(Y))
  fold_list <- split(ind, cut(seq_along(ind), breaks = nfold, labels = FALSE))# Split data according to fold ID
  est.cross <- matrix(ncol = 7, nrow = 0)

  A <- scale(M) %*% t(scale(M))/ncol(M)
  eig_A <- eigen(A, symmetric = TRUE)
  U <- eig_A$vectors
  W <- sqrt(nrow(M)) * U[, 1:3]

  for (i in 1:nfold) {
    hold_out_ind <- fold_list[[i]]
    filter_result <- filter_type2_nonmediator(X[-hold_out_ind], M[-hold_out_ind, ], Y[-hold_out_ind], K, covariates_demo[-hold_out_ind, ], W[-hold_out_ind, ])
    result_i <- estimation(X[hold_out_ind], M[hold_out_ind, ], Y[hold_out_ind], K, covariates_demo[hold_out_ind, ], W[hold_out_ind, ], filter_result)
    est.cross <- rbind(est.cross, result_i)
  }

  result <- apply(est.cross, 2, mean)
  time.total <- as.numeric(Sys.time() - start.time, units = "secs")
  result <- c(result, 'time' = time.total)
  return(result)
}



est_procedure=function(x, M, y, K=0.01, adjustGRM = TRUE, Reg.model = 'ols'){
  p <- ncol(M)
  N <- length(y)
  P = mean(y)
  w = ifelse(y==1, 1, P*(1-K)/K/(1-P)) # weight in IPW

  getAlpha <- function(Var){
    model <- lm(Var ~ x, weights = w)
    alpha.est <- summary(model)$coefficients[2,'Estimate']
    p_value <- summary(model)$coefficients[2,'Pr(>|t|)']
    return(list(alpha.est, p_value))
  }
  temp <- apply(M, 2, getAlpha)
  temp2 <- sapply(temp, function(x) sapply(x, function(y) y[[1]]))
  alpha.est <- temp2[1, ]
  p_value_alpha <- temp2[2, ]

  bias = p*sum((w*x)^2)/(sum(w*(x)^2))^2
  alpha2.est = sum(alpha.est^2) - bias # estimate of alpha^T alpha after bias correction

  sig_alpha <- est_sig_BH(p_value_alpha)
  sig2_a.est <- alpha2.est/p/mean(sig_alpha)

  # estimate residual and alpha
  A <- scale(M) %*% t(scale(M)) / p
  eig_A <- eigen(A, symmetric = TRUE)
  U <- eig_A$vectors
  W <- sqrt(N) * U[, 1:3] # top 3 principal components, check if M is standardized or not here!!!!!!! seems like in simulation, the effect of not standardizing M is not significant.

  residual.est = apply(M, 2, function(Var) {residuals(lm(Var~x + W, weights = w))})
  residual.est = residual.est[, sig_alpha]
  #residual.est = scale(residual.est)
  # estimate sig2_11
  covariates <- cbind(1, cbind(x, W))
  pcgc.est = pcgc(residual.est, y, covariates=covariates, P, K, NA, adjustGRM, Reg.model)
  sig2_11.est = pcgc.est$sig2_g

  # fixed effect glm
  probit.fit = glm(y ~ x + W, family = binomial(link='probit'), weights = w)
  r.est = probit.fit$coefficients[2]

  # estimate of Rmed and total variance
  varl.est = r.est^2 + sig2_a.est*sig2_11.est + 1
  Rmed.est = sig2_a.est*sig2_11.est/varl.est
  V.est = Rmed.est/(Rmed.est + r.est^2/varl.est)
  return(c('Rmed.est'=unname(Rmed.est), 'V.est'=unname(V.est), 'varl.est'= unname(varl.est),
           'sig2_11.est'=unname(sig2_11.est), 'sig2_a.est'=unname(sig2_a.est), 'gamma.est'=unname(r.est)))
}
