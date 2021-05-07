#' Generate a data frame with 10 covariates following a specified distribution
#'
#' Our simulations involve 10 covariates, which are either 1) 10 independent
#' standard normal covariates, 2) 10 standard normal covariates that have pairwise
#' correlations of 0.25, 3) 5 binary Bernouilli covariates with success probability
#' 0.5 and 5 independent standard normal covariates, 4) 9 binary Bernouilli
#' covariates with success probability 0.5 and 1 standard normal covariate, or
#' 5) 10 binary Bernouilli covariates with success probability 0.5. These are
#' abbreviated by "Ind Norm", "Corr Norm", "5 Bern", "9 Bern", and "10 Bern",
#' respectively, and can also be referred to by their number (1-5).
#'
#' @inheritParams calcAlpha0
#' @param n The number of rows or observations to be generated
#'
#' @return Returns a data frame with 10 covariates of the specified distributions
#' and n rows.
#'
#' @export
#' @importFrom stats rnorm rbinom

generateCovariates <- function(cov_dist, n = 5000) {
  # Error checks
  if (!is.wholenumber(n) || length(n) > 1 || n < 1) {
    stop("'n' must be an integer greater or equal to 1.")
  }
  if ( !(cov_dist %in% 1:5 ||
         cov_dist %in% c('Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'))) {
    stop("'cov_dist' must be an integer from 1 to 5 or one of the following: \n
         'Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'.")
  }

  # Create an empty data frame
  data <- data.frame(matrix(NA, nrow = n, ncol = 10,
                            dimnames = list(NULL,
                                            c("X1", "X2", "X3", "X4", "X5", "X6",
                                              "X7", "X8", "X9", "X10"))))

  # Fill in the covariates based on the specified distribution
  if (cov_dist == "Ind Norm" | cov_dist == 1) {
    data[, 1:10] <- matrix(rnorm(n * 10), ncol = 10)
  } else if (cov_dist == "Corr Norm" | cov_dist == 2) {
    Sigma <- matrix(.25, nrow = 10, ncol = 10)
    diag(Sigma) <- 1
    data[, 1:10] <- MASS::mvrnorm(n, rep(0, 10), Sigma)
  } else if (cov_dist == "5 Bern" | cov_dist == 3) {
    data[, 1:5] <- matrix(rbinom(5 * n, 1, .5), ncol = 5)
    data[, 6:10] <- matrix(rnorm(5 * n), ncol = 10)
  } else if (cov_dist == "9 Bern" | cov_dist == 4) {
    data[, 1:9] <- matrix(rbinom(9 * n, 1, .5), ncol = 9)
    data[, 10] <- rnorm(n)
  } else {
    data[, 1:10] <- matrix(rbinom(10 * n, 1, .5), ncol = 10)
  }

  return(data)
}
