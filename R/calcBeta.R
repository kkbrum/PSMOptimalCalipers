#' Calculate the needed value of beta to achieve a specified risk difference
#'
#' For the specified outcome regression, a value of beta must be specified as the
#' log-odds ratio relating the treatment to the outcome. This beta will have a
#' consequence on the true risk difference and thus beta must be calculated accordingly.
#' It is not analytically tractable and thus we estimate it using simulation and the
#' bisection method.
#'
#' @inheritParams calcAlpha0
#' @param true_risk_diff The desired treated risk minus control risk
#' @param alpha_0_outcome The calculated value alpha_0_outcome, obtained from calcAlpha0()
#'
#' @return The value of beta that corresponds to the desired risk difference for
#' the given covariate distribution.
#'
#' @export

calcBeta <- function(true_risk_diff, cov_dist, alpha_0_outcome, tol = 0.0001) {
  # Error checks
  if (!is.numeric(true_risk_diff) || length(true_risk_diff) > 1) {
    stop("'true_risk_diff' must be a numeric.")
  }
  if ( !(cov_dist %in% 1:5 ||
         cov_dist %in% c('Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'))) {
    stop("'cov_dist' must be an integer from 1 to 5 or one of the following: \n
         'Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'.")
  }
  if (!is.numeric(alpha_0_outcome) || length(alpha_0_outcome) > 1) {
    stop("'alpha_0_outcome' must be a numeric.")
  }
  if (!is.numeric(tol) || length(tol) > 1 || tol <= 0) {
    stop("'tol' must be a number greater than 0.")
  }

  # Initialization
  upper <- 0
  lower <- -2
  beta <- mean(c(lower, upper))

  # How many datasets to create in each iteration
  N <- 10
  risk_diff <- rep(NA, N)
  # Run bisection method
  while ((upper - lower) > tol) {
    for (i in 1:N) {
      # Size of each dataset
      n <- 10000
      # Generate covariates
      data <- generateCovariates(cov_dist, n)
      # Given constants
      alpha_L <- log(1.1)
      alpha_M <- log(1.25)
      alpha_H <- log(1.5)
      alpha_VH <- log(2)
      # For each subject, compute both the probability of the outcome if treated and if untreated
      data$p_outcome_1 <- 1 / (1 + exp(- (alpha_0_outcome + beta +
                                            alpha_L * (data$X1 + data$X2 + data$X3) +
                                            alpha_M * (data$X4 + data$X5 + data$X6) +
                                            alpha_H * (data$X7 + data$X8 + data$X9) +
                                            alpha_VH * data$X10)))
      data$p_outcome_0 <- 1 / (1 + exp(- (alpha_0_outcome +
                                            alpha_L * (data$X1 + data$X2 + data$X3) +
                                            alpha_M * (data$X4 + data$X5 + data$X6) +
                                            alpha_H * (data$X7 + data$X8 + data$X9) +
                                            alpha_VH * data$X10)))
      # Record average difference in risks
      risk_diff[i] <- mean(data$p_outcome_1) - mean(data$p_outcome_0)
    }
    est_risk_diff <- mean(risk_diff)
    # Shrink interval of values for beta based on results
    if (est_risk_diff > true_risk_diff) {
      upper <- beta
    } else {
      lower <- beta
    }
    # Take the bisection as new estimate
    beta <- mean(c(lower, upper))
  }
  return(beta)
}
