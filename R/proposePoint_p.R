#' Proposes a new point in the M_p-extension
#' @param x Current point
#' @param B Matrix of sampling directions
#' @param p Relaxation limit
#'
#' @return Returns a new proposed vector in the M_p-extension
#' @export
proposePoint_p <- function(x, B, p){
  b_idx = sample(1:ncol(B), 1, replace = TRUE)
  b = B[, b_idx]
  amax = floor(min(((x+p)/abs(b))[which(b<0)]))
  amin = -floor(min(((x+p)/abs(b))[which(b>0)]))

  xmin = x + amin*b
  a = sample(0:(amax-amin), 1, replace = TRUE)
  x.star = xmin + a*b

  return(x.star)
}
