#' Calculate the needed value of alpha_0_outcome or alpha_0_treat
#'
#' The form of the treatment and outcome regressions require values of
#' alpha_0_outcome or alpha_0_treat. In order to achieve a desired outcome
#' or treatment rate with a given covariate distribution, this value must
#' be calculated. It is not analytically solvable and thus it is approximated
#' using a simulation and bisection method.
#'
#' @param desired_prop The desired event rate
#' @param cov_dist The covariate distribution. An integer 1-5 or one of the following:
#'   'Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'
#' @param tol The maximum difference between the upper and lower estimate
#'
#' @return A numeric to be used as alpha_0_treat or alpha_0_outcome.
#'
#' @export

calcAlpha0 <- function(desired_prop, cov_dist, tol = 0.0001) {
  # Error checks
  if (!is.numeric(desired_prop) || length(desired_prop) != 1 ||
      desired_prop < 0 || desired_prop > 1) {
    stop("'desired_prop' must be a number between 0 and 1.")
  }
  if ( !(cov_dist %in% 1:5 ||
         cov_dist %in% c('Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'))) {
    stop("'cov_dist' must be an integer from 1 to 5 or one of the following: \n
         'Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'.")
  }
  if (!is.numeric(tol) || length(tol) > 1 || tol <= 0) {
    stop("'tol' must be a number greater than 0.")
  }

  # This would give a treated or outcome proportion of 1
  upper <- 0
  # This should be -Inf but by -100 we have a proportion of essentially 0
  lower <- -100
  # Start with a reasonable estimate
  alpha_0 <- log(desired_prop / (1 - desired_prop))

  # How many datasets to create in each iteration
  N <- 10
  # Size of each dataset
  n <- 10000
  prop <- rep(NA, N)
  # Run bisection method
  while ((upper - lower) > tol) {
    for (i in 1:N) {
      # Generate covariates
      data <- generateCovariates(cov_dist, n)
      # Calculate probability of treatment or outcome if untreated
      alpha_L <- log(1.1)
      alpha_M <- log(1.25)
      alpha_H <- log(1.5)
      alpha_VH <- log(2)
      data$p_0 <- 1 / (1 + exp(- (alpha_0 +
        alpha_L * (data$X1 + data$X2 + data$X3) +
        alpha_M * (data$X4 + data$X5 + data$X6) +
        alpha_H * (data$X7 + data$X8 + data$X9) +
        alpha_VH * data$X10)))
      # Record mean baseline risk
      prop[i] <- mean(data$p_0)
    }
    current_prop <- mean(prop)
    # Shrink interval of values for alpha_0 based on results
    if (current_prop > desired_prop) {
      upper <- alpha_0
    } else {
      lower <- alpha_0
    }
    # Take the bisection as new estimate
    alpha_0 <- mean(c(lower, upper))
  }

  return(alpha_0)
}
