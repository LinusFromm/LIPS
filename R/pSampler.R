#' p-samplers
#' @param A Configuration Matrix
#' @param y (Corrupted/aggregated) Observation Vector
#' @param B Set of moves used (p-Markov, CPLB or Full Markov)
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#' @param p Relaxation limit
#' @param loggamma Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param xStart Starting points for each chain. Needs to be fed in as a matrix with number of rows = number of chains
#' @param nSample Number of samples
#' @param nBurnin Number of burnin steps
#' @param chainID Which number chain?
#' @param thinning Thinning parameter
#'
#' @export
pSampler <- function(A,
                     y,
                     B,
                     Model = "uniform",
                     lambda = NULL,
                     p = 0,
                     loggamma = 0,
                     w = 0,
                     xStart = NULL,
                     nSample = 1e+05,
                     nBurnin = 1e+04,
                     chainID = 1,
                     thinning = 1){
  r = nrow(A)
  c = ncol(A)

  if(length(p) == 1){
    p = rep(p, c)
  }

  if(round(nSample/thinning) != nSample/thinning){
    stop("thinning paramter must divide nSample parameter!")
  }

  x = matrix(NA, ncol = c+3, nrow = nSample/thinning)
  x[, c+3] = (chainID-1)*nSample + 1:(nSample/thinning)
  x[, c+2] = 1:(nSample/thinning)
  x[, c+1] = rep(chainID, nSample/thinning)

  if(is.null(xStart)){
    xStart = lpSolve::lp(direction = "min",
                          objective.in = sample(0:1, c, replace = TRUE),
                          const.mat = A,
                          const.rhs = y,
                          const.dir = "=",
                          all.int = TRUE)$solution
  }

  xCurrent = xStart

  moveIndices = sample(1:ncol(B), nBurnin, replace = TRUE)
  for(iiii in 1:nBurnin){
    moveIdx = moveIndices[iiii]
    xProposal = proposePoint(xCurrent, moveIdx, B, extension = "p", p)
    alpha = acceptanceExt(xCurrent, xProposal, loggamma, w, Model, lambda)

    if(stats::runif(1) < exp(alpha)){
      xCurrent = xProposal
    }
  }

  moveIndices = sample(1:ncol(B), nSample, replace = TRUE)
  for(iiiii in 1:nSample){
    moveIdx = moveIndices[iiiii]
    xProposal = proposePoint(xCurrent, moveIdx, B, extension = "p", p)
    alpha = acceptanceExt(xCurrent, xProposal, loggamma, w, Model, lambda)

    if(stats::runif(1) < exp(alpha)){
      xCurrent = xProposal
    }

    if(iiiii%%thinning == 0){
      x[iiiii/thinning, 1:c] = xCurrent
    }
  }
  return(x = x)
}
