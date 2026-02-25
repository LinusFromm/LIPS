#' Computes a column partition lattice basis
#' @param A Configuration matrix
#' @param CPLB_idx Index of columns of A_2
#'
#' @return Returns a column partition lattice basis as matrix
#' @export
computeCPLB <- function(A, CPLB_idx){
  CPLB_idx <- as.integer(CPLB_idx)

  A_1 = A[,-CPLB_idx]
  A_2 = A[,CPLB_idx]

  if(abs(det(A_1)) != 1){
    stop("Column partition lattice basis can only be calculated for unimodular matrix A_1")
  }

  LB_1 = rbind(-solve(A_1)%*%A_2)
  LB_2 = diag(length(CPLB_idx))

  LB = matrix(0, nrow = ncol(A), ncol = length(CPLB_idx))
  LB[CPLB_idx,] = LB_2
  LB[-CPLB_idx,] = LB_1

  return(round(LB))
}
