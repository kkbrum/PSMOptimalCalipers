#' Conduct propensity score matching with a caliper and calculate the treatment effect
#'
#' The logit of the propensity scores are used to match control units to treated
#' units that are within the specified caliper. The risk difference (risk in treated
#' group minus risk in control group) is then calculated on the results matched sample.
#'
#' @param data A data frame containing at least the columns:
#'   'X1', ..., 'X10', 'Z'
#' @param gamma The caliper width in terms of pooled standard deviations
#'
#' @return The risk difference of the matched treated and control units.
#'
#' @export
#' @importFrom stats glm binomial predict var

calcMatchedRiskDiff <- function(data, gamma) {
  # Error checks
  if (! all(c("Z", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10") %in%
            names(data))) {
    stop('"data" must have at least columns with the following names: \n
    "Z", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"')
  }
  if (! is.numeric(gamma) || length(gamma) > 1 || gamma < 0) {
    stop("'gamma' must be a number greater or equal to 0.")
  }

  # Calculate the logit of the propensity scores
  logit_pscore_model <- glm(Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                            data = data, family = binomial(link = "logit"))
  data$logit_pscore <- predict(logit_pscore_model)
  # Convert the gamma value to a caliper
  caliper <- gamma * sqrt((var(data$logit_pscore[data$Z == 0]) +
                             var(data$logit_pscore[data$Z == 1])) / 2)
  # Run propensity score matching using the logit of the propensity scores
  # The caliper here is assumed to be in standard deviations of the logit
  # propensity scores across the sample, but we need the pooled version, and thus
  # we correct for that
  match <- Matching::Match(Tr = data$Z, X = matrix(data$logit_pscore, ncol = 1),
                           caliper = caliper / sqrt(var(data$logit_pscore)),
                           replace = FALSE, ties = FALSE, version = "fast")
  # Calculate the risk difference of the matched sample
  risk_diff <- mean(data$binary_outcome[match$index.treated]) -
    mean(data$binary_outcome[match$index.control])
  return(risk_diff)
}
