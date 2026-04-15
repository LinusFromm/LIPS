#' Calculate acceptance probability
#' @param xCurrent Current state of MCMC
#' @param xProposal Proposed new state of the MCMC
#' @param loggamma Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#'
acceptanceExt <- function(xCurrent, xProposal, loggamma, w, Model, lambda = NULL){
  xCurrent <- as.numeric(xCurrent)
  xProposal <- as.numeric(xProposal)
  loggamma <- as.numeric(loggamma)
  w <- as.numeric(w)
  lambda <- as.numeric(lambda)

  if(Model == "uniform"){
    if(isOutside(xCurrent) && isOutside(xProposal)){
      alpha = min(0, - w * (sum(abs(xProposal[which(xProposal < 0)]))-sum(abs(xCurrent[which(xCurrent < 0)]))))
    } else if(!isOutside(xCurrent) && isOutside(xProposal)){
      alpha = min(0, loggamma - w*sum(abs(xProposal[which(xProposal < 0)])))
    } else if(isOutside(xCurrent) && !isOutside(xProposal)){
      alpha = min(0, w*sum(abs(xCurrent[which(xCurrent < 0)]))-loggamma)
    } else {
      alpha = 0
    }
  } else if(Model == "poisson"){
    if(isOutside(xCurrent) && isOutside(xProposal)){
      p.current = -w*sum(abs(xCurrent[which(xCurrent < 0)])) + sum(stats::dpois(xCurrent[which(xCurrent >= 0)],
                                                                                 lambda = lambda[which(xCurrent >= 0)],
                                                                                 log = TRUE))
      p.proposal = -w*sum(abs(xProposal[which(xProposal < 0)])) + sum(stats::dpois(xProposal[which(xProposal >= 0)],
                                                                                        lambda = lambda[which(xProposal >= 0)],
                                                                                        log = TRUE))
      alpha = min(0, p.proposal - p.current)
    } else if(!isOutside(xCurrent) && isOutside(xProposal)){
      p.current = sum(stats::dpois(xCurrent, lambda = lambda, log = TRUE))
      p.proposal =  -w*sum(abs(xProposal[which(xProposal < 0)])) + sum(stats::dpois(xProposal[which(xProposal >= 0)],
                                                                                     lambda = lambda[which(xProposal >= 0)],
                                                                                     log = TRUE))

      alpha = min(0, loggamma + p.proposal - p.current)
    } else if(isOutside(xCurrent) && !isOutside(xProposal)){
      p.current = -w * sum(abs(xCurrent[which(xCurrent < 0)])) + sum(stats::dpois(xCurrent[which(xCurrent >= 0)],
                                                                                   lambda = lambda[which(xCurrent >= 0)],
                                                                                   log = TRUE))
      p.proposal = sum(stats::dpois(xProposal, lambda = lambda, log = TRUE))

      alpha = min(0, p.proposal - p.current - loggamma)
    } else {
      p.current = sum(stats::dpois(xCurrent, lambda = lambda, log = TRUE))
      p.proposal = sum(stats::dpois(xProposal, lambda = lambda, log = TRUE))

      alpha = min(0, p.proposal - p.current)
    }
  }
  return(alpha)
}
