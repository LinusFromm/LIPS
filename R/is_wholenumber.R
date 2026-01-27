#' Checking for integer values
#'
#' This function tests whether all entries in a vector are integers (up to rounding error)
#' @param x A numerical vector
#' @param tol Tolerance epsilon for being away from a whole number
#' @export
#' @return Returns TRUE if all entries of x are whole numbers; FALSE otherwise
#' @examples
#' is_wholenumber(c(1,2,4.1))

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)  all(abs(x - round(x)) < tol)
