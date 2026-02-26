#' Propose a new point in the hypercube extension
#' @param x Current point
#' @param B Sampling directions (Should be a column partition lattice basis)
#' @param CPLB_idx Index of the A_2 matrix
#' @param a Bounds of the hypercuboid
#'
#' @return Returns a new proposed point in the hypercube extension
#' @export
proposePoint_hypercube <- function(x, B, CPLB_idx, a){
  idx = sample(1:ncol(B), 1, replace = TRUE)
  b = B[,idx]
  xmin = x - x[CPLB_idx[idx]]*b

  c = sample(0:a[idx], 1, replace = TRUE)
  return(xmin + c*b)
}
