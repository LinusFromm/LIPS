#' Calculate acceptance probability
#' @param x.current Current state of MCMC
#' @param x.proposal Proposed new state of the MCMC
#' @param ldelta Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#'
acceptanceExt <- function(x.current, x.proposal, ldelta, w, Model, lambda = NULL){
  x.current <- as.numeric(x.current)
  x.proposal <- as.numeric(x.proposal)
  ldelta <- as.numeric(ldelta)
  w <- as.numeric(w)
  lambda <- as.numeric(lambda)

  if(Model == "uniform"){
    if(isOutside(x.current) && isOutside(x.proposal)){
      alpha = min(0, - w * (sum(abs(x.proposal[which(x.proposal < 0)]))-sum(abs(x.current[which(x.current < 0)]))))
    } else if(!isOutside(x.current) && isOutside(x.proposal)){
      alpha = min(0, ldelta - w*sum(abs(x.proposal[which(x.proposal < 0)])))
    } else if(isOutside(x.current) && !isOutside(x.proposal)){
      alpha = min(0, w*sum(abs(x.current[which(x.current < 0)]))-ldelta)
    } else {
      alpha = 0
    }
  } else if(Model == "poisson"){
    if(isOutside(x.current) && isOutside(x.proposal)){
      p.current = -w*sum(abs(x.current[which(x.current < 0)])) + sum(stats::dpois(x.current[which(x.current >= 0)],
                                                                                 lambda = lambda[which(x.current >= 0)],
                                                                                 log = TRUE))
      p.proposal = -w*sum(abs(x.proposal[which(x.proposal < 0)])) + sum(stats::dpois(x.proposal[which(x.proposal >= 0)],
                                                                                        lambda = lambda[which(x.proposal >= 0)],
                                                                                        log = TRUE))
      alpha = min(0, p.proposal - p.current)
    } else if(!isOutside(x.current) && isOutside(x.proposal)){
      p.current = sum(stats::dpois(x.current, lambda = lambda, log = TRUE))
      p.proposal =  -w*sum(abs(x.proposal[which(x.proposal < 0)])) + sum(stats::dpois(x.proposal[which(x.proposal >= 0)],
                                                                                     lambda = lambda[which(x.proposal >= 0)],
                                                                                     log = TRUE))

      alpha = min(0, ldelta + p.proposal - p.current)
    } else if(isOutside(x.current) && !isOutside(x.proposal)){
      p.current = -w * sum(abs(x.current[which(x.current < 0)])) + sum(stats::dpois(x.current[which(x.current >= 0)],
                                                                                   lambda = lambda[which(x.current >= 0)],
                                                                                   log = TRUE))
      p.proposal = sum(stats::dpois(x.proposal, lambda = lambda, log = TRUE))

      alpha = min(0, p.proposal - p.current - ldelta)
    } else {
      p.current = sum(stats::dpois(x.current, lambda = lambda, log = TRUE))
      p.proposal = sum(stats::dpois(x.proposal, lambda = lambda, log = TRUE))

      alpha = min(0, p.proposal - p.current)
    }
  }
  return(alpha)
}
