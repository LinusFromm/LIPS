#' Propose a new point in the hypercube extension (Gibbs)
#' @param x Current point
#' @param B Sampling directions (Should be a column partition lattice basis)'
#' @param ldelta log(delta) as by definition of this sampler and extension
#' @param w Weight in exponential as described in definition of extension distribution
#' @param CPLB_idx Index of the A_2 matrix
#' @param a Bounds of the hypercuboid
#'
#' @return Returns a new proposed point in the hypercube extension
#' @export
proposePoint_hypercube_gibbs <- function(x, B, ldelta, w, CPLB_idx, a){
  idx = sample(1:ncol(B), 1, replace = TRUE)
  b = B[,idx]
  xmin = x - x[CPLB_idx[idx]]*b
  update_idx = which(b != 0)

  x.matrix = round(t(mapply(seq, from = xmin, by = b, length.out = a[idx]+1)))
  p_in = as.integer(apply(x.matrix, 2, function(col) all(col >= 0)))
  p_out = apply(x.matrix, 2, function(col) ifelse(any(col < 0), exp(-w*sum(abs(col[which(col < 0)]))), 0))

  R_in = sum(p_in)
  R_out = sum(p_out)

  p = (p_in+p_out)/apply(x.matrix, 2, function(col) ifelse(any(col < 0), R_in*exp(-ldelta) + R_out, R_in + exp(ldelta)*R_out))

  c = sample(0:a[idx], 1, replace = TRUE, prob = p/max(p))
  return(xmin + c*b)
}
