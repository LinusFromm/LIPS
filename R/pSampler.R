#' p-samplers
#' @param A Configuration Matrix
#' @param y (Corrupted/aggregated) Observation Vector
#' @param B Set of moves used (p-Markov, CPLB or Full Markov)
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#' @param p Relaxation limit
#' @param ldelta Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param x.start Starting points for each chain. Needs to be fed in as a matrix with number of rows = number of chains
#' @param n.sample Number of samples
#' @param n.burnin Number of burnin steps
#' @param chain.id Which number chain?
#' @param thinning Thinning parameter
#'
#' @export
pSampler <- function(A,
                     y,
                     B,
                     Model = "uniform",
                     lambda = NULL,
                     p = 0,
                     ldelta = 0,
                     w = 0,
                     x.start = NULL,
                     n.sample = 1e+05,
                     n.burnin = 1e+04,
                     chain.id = 1,
                     thinning = 1){
  r = nrow(A)
  c = ncol(A)

  if(length(p) == 1){
    p = rep(p, c)
  }

  if(round(n.sample/thinning) != n.sample/thinning){
    stop("thinning paramter must divide n.sample parameter!")
  }

  x = matrix(NA, ncol = c+3, nrow = n.sample/thinning)
  x[, c+3] = (chain.id-1)*n.sample + 1:(n.sample/thinning)
  x[, c+2] = 1:(n.sample/thinning)
  x[, c+1] = rep(chain.id, n.sample/thinning)

  if(is.null(x.start)){
    x.start = lpSolve::lp(direction = "min",
                          objective.in = sample(0:1, c, replace = TRUE),
                          const.mat = A,
                          const.rhs = y,
                          const.dir = "=",
                          all.int = TRUE)$solution
  }

  x.current = x.start

  moveIndices = sample(1:ncol(B), n.burnin, replace = TRUE)
  for(iiii in 1:n.burnin){
    moveIdx = moveIndices[iiii]
    x.proposal = proposePoint(x.current, moveIdx, B, extension = "p", p)
    alpha = acceptanceExt(x.current, x.proposal, ldelta, w, Model, lambda)

    if(stats::runif(1) < exp(alpha)){
      x.current = x.proposal
    }
  }

  moveIndices = sample(1:ncol(B), n.sample, replace = TRUE)
  for(iiiii in 1:n.sample){
    moveIdx = moveIndices[iiiii]
    x.proposal = proposePoint(x.current, moveIdx, B, extension = "p", p)
    alpha = acceptanceExt(x.current, x.proposal, ldelta, w, Model, lambda)

    if(stats::runif(1) < exp(alpha)){
      x.current = x.proposal
    }

    if(iiiii%%thinning == 0){
      x[iiiii/thinning, 1:c] = x.current
    }
  }
  return(x = x)
}
