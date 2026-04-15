#' Computes a column partition lattice basis
#' @param A Configuration matrix
#' @param A2Idx Index of columns of A2
#'
#' @return Returns a column partition lattice basis as matrix
#' @export
computeCPLB <- function(A, A2Idx){
  A2Idx <- as.integer(A2Idx)

  A1 = A[,-A2Idx]
  A2 = A[,A2Idx]

  if(abs(det(A1)) != 1){
    stop("Column partition lattice basis can only be calculated for unimodular matrix A1")
  }

  LB1 = rbind(-solve(A1)%*%A2)
  LB2 = diag(length(A2Idx))

  LB = matrix(0, nrow = ncol(A), ncol = length(A2Idx))
  LB[A2Idx,] = LB2
  LB[-A2Idx,] = LB1

  return(round(LB))
}
