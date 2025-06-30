#' @title Cross-Fitted R-squared Based Mediation Effect Estimation
#' @description This function estimates the mediation effect using cross-fitting.
#'
#' @param X Numeric vector. Exposure variable.
#' @param M Numeric matrix. Mediators (rows = samples, columns = mediators).
#' @param Y Numeric vector. Binary outcome variable (0/1).
#' @param K Numeric. Proportion of cases in the population.
#' @param covariates Numeric matrix. Covariates included in the model.
#' @param nfold Integer. Number of cross-fitting folds. Default is 2.
#' @param fdr_method Str. FDR control method, can be either 'storey' or any options provided in function p.adjust(). Default is 'BH'.
#' @param fdr_cutoff Numeric. FDR cutoff value. Default is 0.01.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A named list containing:
#' \describe{
#'   \item{Rmed.est}{Estimated mediation effect.}
#'   \item{Q.est}{Estimated Q measure.}
#'   \item{varl.est}{Total variance estimate.}
#'   \item{sig2_11.est}{Estimated variance component for indirect effects.}
#'   \item{sig2_a.est}{Estimated variance component for exposure effects.}
#'   \item{gamma.est}{Estimated direct effect of exposure on outcome.}
#'   \item{sig_number}{Number of selected significant mediators.}
#'   \item{time}{Total computation time in seconds.}
#' }
#'
#' @examples
#' data <- data_generate(p = 2000, N = 500, K = 0.05, r = 0.5, res.var = 0.2,
#'                       pi_11 = 0.09, pi_alpha = 0.3, pi_beta = 0.3,
#'                       sig2_a = 0.3, sig2_11 = 0.5, sig2_01 = 0.3, seed = 23)
#' X <- data$x
#' M <- data$M
#' Y <- data$y
#' K <- 0.05
#'
#' result <- r2_estimation_cf(X, M, Y, K)
#' print(result)
#'
#' @export

# n-fold cross fitting
r2_estimation_cf <- function(X, M, Y, K, covariates = NULL, nfold = 2, fdr_method = 'BH', fdr_cutoff = 0.01, seed = 123) {
  start.time <- Sys.time()
  set.seed(seed)
  ind <- sample(length(Y))
  fold_list <- split(ind, cut(seq_along(ind), breaks = nfold, labels = FALSE)) # Split data according to fold ID
  est.cross <- matrix(ncol = 7, nrow = 0)

  A <- scale(M) %*% t(scale(M))/ncol(M)
  eig_A <- eigen(A, symmetric = TRUE)
  U <- eig_A$vectors
  W <- sqrt(nrow(M)) * U[, 1:3]

  for (i in 1:nfold) {
    hold_out_ind <- fold_list[[i]]
    covariates_filter <- if (!is.null(covariates)) covariates[-hold_out_ind, ] else NULL
    covariates_est <- if (!is.null(covariates)) covariates[hold_out_ind, ] else NULL

    filter_result <- filter_type2_nonmediator(X[-hold_out_ind], M[-hold_out_ind, ], Y[-hold_out_ind], K, covariates_filter, W[-hold_out_ind, ], fdr_method, fdr_cutoff)
    result_i <- estimation(X[hold_out_ind], M[hold_out_ind, ], Y[hold_out_ind], K, covariates_est, W[hold_out_ind, ], filter_result)
    est.cross <- rbind(est.cross, result_i)
  }

  result <- apply(est.cross, 2, mean)
  time.total <- as.numeric(Sys.time() - start.time, units = "secs")
  result <- c(result, 'time' = time.total)
  return(as.list(result))
}


# pcgc algorithm
pcgc = function(genotype.stand, phenotype, covariates=NA, P, K, P_cond=NA, adjustGRM = TRUE){
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
  model=lm(Zij~-1+Gij)
  sig2_g=coefficients(model)

  return(list('sig2_g' = sig2_g, 'sig2_t' = sig2_t, 'P_cond' = P_cond, 'gamma' = gamma))
}

# step 1: mediatior selection and estimation of sig_a
filter_type2_nonmediator <- function(X, M, Y, K, covariates_demo = NULL, W, fdr_method = 'BH', fdr_cutoff = 0.01){
  P <- mean(Y)
  N <- length(Y)
  p <- ncol(M)
  w <- ifelse(Y == 1, 1, P * (1 - K) / K / (1 - P))
  if (is.null(covariates_demo)) {
    covariates <- W
  } else {
    covariates <- cbind(covariates_demo, W)
  }

  residual.est <- apply(M, 2, function(Var) {residuals(lm(Var ~ X + covariates, weights = w))}) # new residual with mean 0 and sd 1 in population

  ind_cor <- remove_highly_correlated(residual.est, threshold = 0.8)

  getAlpha <- function(Var){
    if (is.null(covariates_demo)) {
      model <- lm(Var ~ X, weights = w)
    } else {
      model <- lm(Var ~ X + covariates_demo, weights = w)
    }
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

  sig_alpha <- est_sig(p_value_alpha, method = fdr_method, cutoff = fdr_cutoff)
  sig2_a.est <- alpha2.est/p/mean(sig_alpha)

  return(list('selected_var' = unname(sig_alpha & ind_cor), 'sig2_a.est' = sig2_a.est))
}

# step 2 and 3: estimation of gamma and sig_11
estimation <- function(X, M, Y, K, covariates_demo = NULL, W, filter_result) {
  P <- mean(Y)
  N <- length(Y)
  p <- ncol(M)
  w <- ifelse(Y == 1, 1, P * (1 - K) / K / (1 - P))
  if (is.null(covariates_demo)) {
    covariates <- cbind(matrix(1, nrow=nrow(M)), W)
  } else {
    covariates <- cbind(covariates_demo, W)
  }

  # fixed effect glm to estimate gamma
  probit.fit <- glm(Y ~ X + covariates, family = binomial(link='probit'), weights = w)
  r.est <- probit.fit$coefficients[2]

  # calculate mediator residual
  residual.est <- apply(M, 2, function(Var) {residuals(lm(Var ~ X + covariates, weights = w))}) # new residual with mean 0 and sd 1 in population

  # get selected mediators and sig2_a from other folds
  sig_alpha <- filter_result$selected_var
  sig2_a.est <- filter_result$sig2_a.est
  if (all(!sig_alpha)) { # no mediation effect if all of the features are type2 non-mediators.
    varl.est <- r.est^2 + 1
    return(c('Rmed.est' = 0, 'Q.est' = 0, 'varl.est' = unname(varl.est),
             'sig2_11.est' = 0, 'sig2_a.est' = unname(sig2_a.est), 'gamma.est' = unname(r.est),
             'sig_number' = sum(sig_alpha)))
  }
  residual.est.filter <- residual.est[, sig_alpha]

  # estimate sig2_11
  pcgc.est <- pcgc(residual.est.filter, Y, covariates = cbind(X, covariates), P = P, K = K, P_cond = NA, adjustGRM=TRUE)
  sig2_11.est <- pcgc.est$sig2_g

  # estimate of Rmed and total variance
  varl.est <- r.est^2 + sig2_a.est*sig2_11.est + 1
  Rmed.est <- sig2_a.est*sig2_11.est/varl.est
  Q.est <- Rmed.est/(Rmed.est + r.est^2/varl.est)

  if (r.est^2 < 0.0001 && sig2_a.est * sig2_11.est < 0.0001) {
    warning("Q_2 can be unstable when total effect size is small.")
  }

  return(c('Rmed.est' = unname(Rmed.est), 'Q.est' = unname(Q.est), 'varl.est' = unname(varl.est),
            'sig2_11.est' = unname(sig2_11.est), 'sig2_a.est' = unname(sig2_a.est), 'gamma.est' = unname(r.est),
            'sig_number' = sum(sig_alpha)))
}

#' @title Jackknife Confidence Interval for Mediation Effect
#' @description Estimates 95\% jackknife confidence intervals for R2-based mediation analysis using cross-fitting.
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#'
#' @param cores Integer. Number of CPU cores to use for parallel computation.
#' @param X Numeric vector. Exposure variable.
#' @param M Numeric matrix. Mediators (rows = samples, columns = mediators).
#' @param Y Numeric vector. Binary outcome variable (0/1).
#' @param K Numeric. Proportion of cases in the population.
#' @param covariates Numeric matrix. Optional covariates to include in the model. Default is NULL.
#' @param nfold Integer. Number of folds for cross-fitting. Default is 2.
#' @param fdr_method Character. Method for false discovery rate control. Default is "BH".
#' @param fdr_cutoff Numeric. FDR cutoff threshold. Default is 0.01.
#'
#' @return A matrix with two rows (2.5\% and 97.5\% quantiles) and columns corresponding to:
#' \itemize{
#'   \item Rmed.est
#'   \item Q.est
#'   \item varl.est
#'   \item sig2_11.est
#'   \item sig2_a.est
#'   \item gamma.est
#'   \item sig_number
#'   \item time
#' }
#'
#' @export

confidence_interval <- function(cores, X, M, Y, K, covariates = NULL, nfold = 2, fdr_method = 'BH', fdr_cutoff = 0.01){
  doParallel::registerDoParallel(cores = cores)
  A <- scale(M) %*% t(scale(M))/ncol(M)
  eig_A <- eigen(A, symmetric = TRUE)
  U <- eig_A$vectors
  W <- sqrt(nrow(M)) * U[, 1:3]
  est.jack <- foreach::foreach(i = 1:length(Y), .combine = rbind) %dopar% {
    tryCatch(
      {
        unlist(r2_estimation_cf(X[-i], M[-i, ], Y[-i], K, covariates[-i, ], nfold, fdr_method, fdr_cutoff))
      },
      error = function(e) {
        message(paste("Error in iteration", i, ":", e$message))
        # Return a row of NA values (or handle in a way consistent with your data)
        return(rep(NA, 8))
      }
    )
  }
  ci.result <- apply(est.jack, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  return(ci.result)
}

#' @title Weak effect plot
#' @description This function visualizes the distribution of p-values for the association between each mediator and a binary outcome. It highlights potential weak mediator-outcome associations when the mediator is also associated with the exposure.
#'
#' @importFrom qvalue qvalue
#' @importFrom ggplot2 ggplot aes geom_histogram geom_hline annotate labs theme_minimal
#'
#' @param X Numeric vector. Exposure variable.
#' @param M Numeric matrix. Mediators (rows = samples, columns = mediators).
#' @param Y Numeric vector. Binary outcome variable (0/1).
#' @param K Numeric. Proportion of cases in the population.
#' @param covariates Numeric matrix. Covariates included in the model.
#'
#' @return A ggplot2 object showing the histogram of p-values for the mediator-outcome associations.
#'
#' @examples
#' data <- data_generate(p = 2000, N = 500, K = 0.05, r = 0.5, res.var = 0.2,
#'                       pi_11 = 0.09, pi_alpha = 0.3, pi_beta = 0.3,
#'                       sig2_a = 0.3, sig2_11 = 0.5, sig2_01 = 0.3, seed = 23)
#' X <- data$x
#' M <- data$M
#' Y <- data$y
#' K <- 0.05
#'
#' result <- r2_estimation_cf(X, M, Y, K)
#' print(result)
#'
#' @export

weak_effect_plot <- function(X, M, Y, K, covariates = NULL){
  P <- mean(Y)
  N <- length(Y)
  if (is.null(covariates)) {
    covariates <- rep(1, N)
  }
  p <- ncol(M)
  w <- ifelse(Y == 1, 1, P * (1 - K) / K / (1 - P))
  A <- scale(M) %*% t(scale(M))/ncol(M)
  eig_A <- eigen(A, symmetric = TRUE)
  U <- eig_A$vectors
  W <- sqrt(nrow(M)) * U[, 1:6]

  # q values for alpha
  getAlpha <- function(Var){
    model <- lm(Var ~ X + covariates, weights = w)
    alpha.est <- summary(model)$coefficients[2,'Estimate']
    p_value <- summary(model)$coefficients[2,'Pr(>|t|)']
    return(list(alpha.est, p_value))
  }

  temp <- apply(M, 2, getAlpha)
  temp2 <- sapply(temp, function(x) sapply(x, function(y) y[[1]]))
  alpha.est <- temp2[1, ]
  p_values_alpha <- temp2[2, ]

  # q values for beta
  getBeta <- function(Var){
    model <- glm(Y ~ Var + X + covariates + W, family = binomial(link='probit'), weights = w)
    beta.est <- summary(model)$coefficients[2,'Estimate']
    p_value <- summary(model)$coefficients[2,'Pr(>|z|)']
    return(list(beta.est, p_value))
  }

  temp <- apply(M, 2, getBeta)
  temp2 <- sapply(temp, function(x) sapply(x, function(y) y[[1]]))
  beta.est <- temp2[1, ]
  p_values_beta <- temp2[2, ]
  data <- data.frame(p_values_beta = p_values_beta, p_values_alpha = p_values_alpha)

  # weak effect plot
  qobj <- qvalue(p_values_beta)
  pi0_estimate <- qobj$pi0  # Estimated proportion of null p-values
  flatten_density <- pi0_estimate * 1

  ggplot(subset(data, p_values_alpha < 0.05), aes(x = p_values_beta)) +
    geom_histogram(aes(y = ..density..), breaks = seq(0, 1, length.out = 21), fill = "skyblue", color = "black", alpha = 0.5) +
    geom_hline(yintercept = flatten_density, linetype = "dashed", color = "red", linewidth = 0.5) +
    annotate("text", x = 0.7, y = flatten_density + 1.5,
             label = paste0("Estimated null density level = ", round(flatten_density, 2)),
             color = "red", size = 4)+
    labs(x = "P-values for Mediator-Outcome Associations",
         y = "Density") +
    theme_minimal()
}

