#' Calculate Acceptance ratio
#' @param xCurrent Current point
#' @param xProposal Proposed point
#' @param loggamma log(delta) used to define epsilon
#' @param w Weight of negativity
#' @param dist Model distribution
#'
#' @return Returns the Acceptance ratio
#' @export
calculateAlpha <- function(xCurrent, xProposal, loggamma, w = 0, dist = "Unif"){

  xCurrent <- as.numeric(xCurrent)
  xProposal <- as.numeric(xProposal)
  loggamma <- as.numeric(loggamma)

  if(dist == "unif"){
    if(isOutside(xCurrent) && isOutside(xProposal)){
      alpha = min(0, w * (sum(abs(xCurrent[which(xCurrent < 0)]))-sum(abs(xProposal[which(xProposal < 0)]))))
    } else if(!isOutside(xCurrent) && isOutside(xProposal)){
      alpha = min(0, loggamma - w*sum(abs(xProposal[which(xProposal < 0)])))
    } else if(isOutside(xCurrent) && !isOutside(xProposal)){
      alpha = min(0, w*sum(abs(xCurrent[which(xCurrent < 0)]))-loggamma)
    } else {
      alpha = 0
    }
  } else if(dist == "no3way"){
    if(isOutside(xCurrent) && isOutside(xProposal)){
      alpha = min(0, w * (sum(abs(xCurrent[which(xCurrent < 0)]))-sum(abs(xProposal[which(xProposal < 0)]))))
    } else if(!isOutside(xCurrent) && isOutside(xProposal)){
      alpha = min(0, loggamma + sum(lfactorial(xCurrent)) - w*sum(abs(xProposal[which(xProposal < 0)])))
    } else if(isOutside(xCurrent) && !isOutside(xProposal)){
      alpha = min(0, w*sum(abs(xCurrent[which(xCurrent < 0)]))-loggamma - sum(lfactorial(xProposal)))
    } else {
      alpha = min(0, -(sum(lfactorial(xProposal)) - sum(lfactorial(xCurrent))))
    }
  }
  return(alpha)
}
