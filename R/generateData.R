#' Generate a data frame containing 10 covariates, a treatment indicator, and a binary outcome
#'
#' This function generates the needed data frame for the simulations,
#' which includes 10 covariates of a specified distribution, a treatment indicator,
#' and a binary outcome variable.
#'
#' @inheritParams calcAlpha0
#' @inheritParams calcBeta
#' @param alpha_0_treat The calculated value of alpha_0_treat, obtained from calcAlpha0()
#' @param beta The calculated value of beta, obtained from calcBeta()
#'
#' @return A data frame with columns for the covariates, treatment indicator, and binary outcome,
#' as well as the logit and regular probabilities of treatment and outcome that were used
#' to generate them.
#'
#' @export
#' @importFrom stats rbinom

generateData <- function(cov_dist, alpha_0_treat, alpha_0_outcome, beta) {
  # Error checks
  if ( !(cov_dist %in% 1:5 ||
         cov_dist %in% c('Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'))) {
    stop("'cov_dist' must be an integer from 1 to 5 or one of the following: \n
         'Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'.")
  }
  if (!is.numeric(alpha_0_treat) || length(alpha_0_treat) > 1) {
    stop("'alpha_0_treat' must be a numeric.")
  }
  if (!is.numeric(alpha_0_outcome) || length(alpha_0_outcome) > 1) {
    stop("'alpha_0_outcome' must be a numeric.")
  }
  if (!is.numeric(beta) || length(beta) > 1) {
    stop("'beta' must be a numeric.")
  }

  n <- 5000
  # Generate covariates
  data <- generateCovariates(cov_dist, n)
  # Generate treatment indicator
  alpha_L <- log(1.1)
  alpha_M <- log(1.25)
  alpha_H <- log(1.5)
  alpha_VH <- log(2)
  data$logit_p_treat <- alpha_0_treat +
    alpha_L * (data$X1 + data$X2 + data$X3) +
    alpha_M * (data$X4 + data$X5 + data$X6) +
    alpha_H * (data$X7 + data$X8 + data$X9) +
    alpha_VH * data$X10
  data$p_treat <- 1 / (1 + exp(- data$logit_p_treat))
  data$Z <- rbinom(n, 1, data$p_treat)
  # Generate binary outcome
  data$logit_p_outcome <- alpha_0_outcome + beta * data$Z +
    alpha_L * (data$X1 + data$X2 + data$X3) +
    alpha_M * (data$X4 + data$X5 + data$X6) +
    alpha_H * (data$X7 + data$X8 + data$X9) +
    alpha_VH * data$X10
  data$p_outcome <- 1 / (1 + exp( - data$logit_p_outcome))
  data$binary_outcome <- rbinom(n, 1, data$p_outcome)
  return(data)
}
