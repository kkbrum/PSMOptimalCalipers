#' Run a simulation with a specified caliper width, risk difference, and covariate distribution
#'
#' This runs a simulation that creates N data sets, each of size 10,000, that have a
#' specified risk difference and covariate distribution. Propensity score matching
#' with a caliper derived from the one given gamma is then conducted and the risk difference
#' in the full and matched data sets is estimated. The average estimated risk difference
#' for the full and matched data sets is reported, as well as the bias reduction of
#' matching and the MSE of the matched estimator. The mean and variance of the time it took
#' to generate 1 data set, match, and calculate risk differences is also reported.
#'
#' @inheritParams calcMatchedRiskDiff
#' @inheritParams calcBeta
#' @inheritParams generateData
#' @param N The number of data sets to generate within the simulation
#'
#' @return A named vector containing the gamma value used, MSE of the matched estimator,
#' the bias reduction of using the matched estimator, the average estimate from the
#' full data, the average estimate from the matched data, the average time to go through
#' one data set, and the variance of the times to go through one data set.
#'
#' @export
#' @importFrom stats var

runSimulation <- function(gamma, N, true_risk_diff, cov_dist,
                          alpha_0_treat, alpha_0_outcome, beta) {
  # Error checks
  if (! is.numeric(gamma) || length(gamma) > 1 || gamma < 0) {
    stop("'gamma' must be a number greater or equal to 0.")
  }
  if (!is.wholenumber(N) || length(N) > 1 || N < 1) {
    stop("'N' must be an integer greater than or equal to 1.")
  }
  if (!is.numeric(true_risk_diff) || length(true_risk_diff) > 1) {
    stop("'true_risk_diff' must be a numeric.")
  }
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

  # Create vectors to store results
  risk_diff <- rep(NA, N)
  crude_risk_diff <- rep(NA, N)
  times <- rep(NA, N)
  # Run the simulation for each data set
  for (i in 1:N) {
    start <- Sys.time()
    data <- generateData(cov_dist = cov_dist, alpha_0_treat = alpha_0_treat,
                         alpha_0_outcome = alpha_0_outcome, beta = beta)
    # Calculate the risk difference in the full data
    crude_risk_diff[i] <- mean(data$binary_outcome[data$Z == 1]) - mean(data$binary_outcome[data$Z == 0])
    # Calculate the risk difference in the matched data set
    risk_diff[i] <- calcMatchedRiskDiff(data, gamma)
    end <- Sys.time()
    times[i] <- end - start
  }
  # Calculate and store various metrics about the N data sets
  MSE <- mean((risk_diff - true_risk_diff)^2)
  bias_crude <- abs(mean(crude_risk_diff) - true_risk_diff)
  bias_ps <- abs(mean(risk_diff) - true_risk_diff)
  mean_crude <- mean(crude_risk_diff)
  mean_ps <- mean(risk_diff)
  bias_reduction <- 100 * (bias_crude - bias_ps) / bias_crude
  time_mean <- mean(times)
  time_var <- var(times)
  return(c(gamma = gamma, MSE = MSE, bias_reduction = bias_reduction,
           mean_crude = mean_crude, mean_ps = mean_ps,
           time_mean = time_mean, time_var = time_var))
}
