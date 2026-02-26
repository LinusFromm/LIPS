#' Calculate Acceptance ratio
#' @param x Current point
#' @param x.star Proposed point
#' @param ldelta log(delta) used to define epsilon
#' @param w Weight of negativity
#' @param dist Model distribution
#'
#' @return Returns the Acceptance ratio
#' @export
calculateAlpha <- function(x, x.star, ldelta, w = 0, dist = "Unif"){

  x <- as.numeric(x)
  x.star <- as.numeric(x.star)
  ldelta <- as.numeric(ldelta)

  if(dist == "unif"){
    if(isOutside(x) && isOutside(x.star)){
      alpha = min(0, w * (sum(abs(x[which(x < 0)]))-sum(abs(x.star[which(x.star < 0)]))))
    } else if(!isOutside(x) && isOutside(x.star)){
      alpha = min(0, ldelta - w*sum(abs(x.star[which(x.star < 0)])))
    } else if(isOutside(x) && !isOutside(x.star)){
      alpha = min(0, w*sum(abs(x[which(x < 0)]))-ldelta)
    } else {
      alpha = 0
    }
  } else if(dist == "no3way"){
    if(isOutside(x) && isOutside(x.star)){
      alpha = min(0, w * (sum(abs(x[which(x < 0)]))-sum(abs(x.star[which(x.star < 0)]))))
    } else if(!isOutside(x) && isOutside(x.star)){
      alpha = min(0, ldelta + sum(lfactorial(x)) - w*sum(abs(x.star[which(x.star < 0)])))
    } else if(isOutside(x) && !isOutside(x.star)){
      alpha = min(0, w*sum(abs(x[which(x < 0)]))-ldelta - sum(lfactorial(x.star)))
    } else {
      alpha = min(0, -(sum(lfactorial(x.star)) - sum(lfactorial(x))))
    }
  }
  return(alpha)
}
