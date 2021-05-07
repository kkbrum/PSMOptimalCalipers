#' Run simulations for various covariate distributions, risk differences, and caliper widths
#'
#' This function runs simulations for a combination of covariate distributions,
#' risk differences and caliper widths and a given treatment proportion and base outcome risk.
#' The function can be parallelized with multiple cores, and reproduced with a seed value.
#'
#' @inheritParams runSimulation
#' @param seed The seed to use for the random number generator
#' @param cores The number of cores to use
#' @param cov_dists A vector of covariate distributions for which to run simulations
#' @param true_risk_diffs A vector of true risk differences for which to run simulations
#' @param gammas A vector of caliper widths in terms of pooled standard deviations for which to run simulations
#' @param p_treat A numeric between 0 and 1 giving the proportion of units that should be treated (constant across simulations)
#' @param base_risk A numeric between 0 and 1 giving the base risk of the outcome if all units untreated
#' @param file Optional file to save workspace to after each set of N simulations
#'
#' @return A list with an entry for each covariate entry. That entry is a list with
#' an entry for each of the true risk differences. This entry is a matrix, with one
#' column for each gamma, which contains the set of simulation results:
#' the gamma value used, MSE of the matched estimator,
#' the bias reduction of using the matched estimator, the average estimate from the
#' full data, the average estimate from the matched data, the average time to go through
#' one data set, and the variance of the times to go through one data set.
#'
#' @examples
#' runSimulations(seed = 64, cores = 1, cov_dists = c("Ind Norm"),
#'   true_risk_diffs = c(-0.05), gammas = c(1), N = 1, p_treat = 0.25,
#'   base_risk = 0.29, file = NULL)
#'
#' @export

runSimulations <- function(seed, cores, cov_dists, true_risk_diffs, gammas, N,
                           p_treat, base_risk, file = NULL) {
  # Error checks
  if (!is.wholenumber(seed) || length(seed) > 1 || seed < 1) {
    stop("'seed' must be an integer greater than or equal to 1.")
  }
  if (!is.wholenumber(cores) || length(cores) > 1 || cores < 1) {
    stop("'cores' must be an integer greater than or equal to 1.")
  }
  if (!is.wholenumber(N) || length(N) > 1 || N < 1) {
    stop("'N' must be an integer greater than or equal to 1.")
  }
  if ( !all(cov_dists %in% 1:5 |
         cov_dists %in% c('Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'))) {
    stop("'cov_dist' must be an integer from 1 to 5 or one of the following: \n
         'Ind Norm', 'Corr Norm', '5 Bern', '9 Bern', '10 Bern'.")
  }
  if (!is.numeric(true_risk_diffs)) {
    stop("'true_risk_diffs' must be a scalar or vector of numerics.")
  }
  if (!is.numeric(gammas) || !all(gammas >= 0)) {
    stop("'gammas' must be a scalar or vector of numerics greater or equal to 0.")
  }
  if (!is.wholenumber(N) || length(N) > 1 || N < 1) {
    stop("'N' must be an integer greater than or equal to 1.")
  }
  if (!is.numeric(p_treat) || length(p_treat) != 1 ||
      p_treat < 0 || p_treat > 1) {
    stop("'p_treat' must be a number between 0 and 1.")
  }
  if (!is.numeric(base_risk) || length(base_risk) != 1 ||
      base_risk < 0 || base_risk > 1) {
    stop("'base_risk' must be a number between 0 and 1.")
  }

  # Set RNG for parallelization
  RNGkind("L'Ecuyer-CMRG")
  # Create cluster for parallelization
  cl <- parallel::makeCluster(cores)
  parallel::clusterSetRNGStream(cl, iseed = seed)

  # Create a list of lists to store the results in
  results <- vector("list", length(cov_dists))
  names(results) <- cov_dists
  for (i in 1:length(cov_dists)) {
    results[[i]] <- vector("list", length(true_risk_diffs))
    names(results[[i]]) <- paste0("RD ", true_risk_diffs)
  }

  # We run simulations for each of the covariate distributions
  for (j in 1:length(cov_dists)) {
    cov_dist <- cov_dists[j]
    alpha_0_treat <- calcAlpha0(p_treat, cov_dist)
    alpha_0_outcome <- calcAlpha0(base_risk, cov_dist)

    # We run simulations for each of the risk differences
    for (i in 1:length(true_risk_diffs)) {
      true_risk_diff <- true_risk_diffs[i]
      beta <- calcBeta(true_risk_diff, cov_dist,
                       alpha_0_outcome)
      results[[j]][[i]] <- parallel::parSapply(cl, X = gammas, runSimulation,
                                     N = N, true_risk_diff = true_risk_diff,
                                     cov_dist = cov_dist, alpha_0_treat = alpha_0_treat,
                                     alpha_0_outcome = alpha_0_outcome, beta = beta)
      if (! is.null(file)) {
        save(results, file = file)
      }
    }
  }

  parallel::stopCluster(cl)

  return(results = results)
}
