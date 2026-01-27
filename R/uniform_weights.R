#' Computes geometric weights for every move in B
#' @param A Configuration matrix
#' @param y Observation vector
#' @param B Matrix containing basis used in the Hit-and-Run sampler. Every column is a move from the basis
#'
#' @return Returns a vector of weights used for the weighted_move_sampler
#' @export
uniform_weights <- function(A,
                            y,
                            B){
  c = ncol(A)
  r = nrow(A)
  m = ncol(B)

  x.max = numeric(c)
  x.min = numeric(c)

  for(i in 1:c){
    x.min[i] = lpSolve::lp(direction = "min", objective.in = diag(c)[i,], const.mat = A, const.dir = "=", const.rhs = y, all.int = TRUE)$solution[i]
    x.max[i] = lpSolve::lp(direction = "max", objective.in = diag(c)[i,], const.mat = A, const.dir = "=", const.rhs = y, all.int = TRUE)$solution[i]
  }

  weights = apply(floor((x.max-x.min)/abs(B)), 2, function(col) mean(col[col < Inf]))
  return(weights)
}
