% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_generation.R
\name{data_generate}
\alias{data_generate}
\title{Generate Synthetic Data for High-Dimensional Mediation Analysis}
\usage{
data_generate(
  p,
  N,
  K,
  r,
  res.var,
  pi_11 = NULL,
  pi_alpha = NULL,
  pi_beta = NULL,
  sig2_a = NULL,
  sig2_11 = NULL,
  sig2_01 = NULL,
  alpha = NULL,
  beta = NULL,
  seed = 123
)
}
\arguments{
\item{p}{Integer. Number of mediators.}

\item{N}{Integer. Total sample size.}

\item{K}{Numeric. Proportion of cases in the population.}

\item{r}{Numeric. Direct effect of exposure on the outcome.}

\item{res.var}{Numeric. Residual variance in the outcome.}

\item{pi_11}{Numeric. Proportion of mediators with a nonzero effect for both alpha and beta.}

\item{pi_alpha}{Numeric. Proportion of mediators with a nonzero alpha effect.}

\item{pi_beta}{Numeric. Proportion of mediators with a nonzero beta effect.}

\item{sig2_a}{Numeric. Variance of alpha effects.}

\item{sig2_11}{Numeric. Variance of beta effects for shared mediators.}

\item{sig2_01}{Numeric. Variance of beta effects for independent mediators.}

\item{alpha}{Numeric vector. Directly specified alpha values (optional).}

\item{beta}{Numeric vector. Directly specified beta values (optional).}

\item{seed}{Integer. Random seed for reproducibility.}
}
\value{
A list containing:
\describe{
  \item{x}{Numeric vector. Exposure values.}
  \item{y}{Binary outcome vector (0/1).}
  \item{liab}{Numeric vector. Underlying liability scores.}
  \item{M}{Matrix. Mediator values.}
  \item{Q}{Q measurement.}
  \item{Rmed}{Estimated mediation effect (theoretical).}
  \item{varl}{Variance of liability.}
  \item{Q.emp}{Empirical Q measurement.}
  \item{Rmed.emp}{Empirical mediation effect.}
  \item{alpha}{Vector of generated alpha values.}
  \item{beta}{Vector of generated beta values.}
}
}
\description{
Simulates exposure, mediators, and binary outcomes with specified variance components. The function generates mediator effects and performs case-control sampling.
}
\examples{
data <- data_generate(p = 1000, N = 100, K = 0.05, r = 0.2, res.var = 0.5,
                      pi_11 = 0.09, pi_alpha = 0.3, pi_beta = 0.3,
                      sig2_a = 0.5, sig2_11 = 0.2, sig2_01 = 0.3, seed = 123)
str(data)

}
