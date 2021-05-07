#' Utility function to check if a number is a whole number
#'
#' @param x The number to check
#' @param tol The tolerance
#'
#' @return Boolean stating whether x is a whole number

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
