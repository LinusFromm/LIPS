#' Extension Samplers
#' @param A Configuration Matrix
#' @param y (Corrupted/aggregated) Observation Vector
#' @param B Set of moves used (p-Markov, CPLB or Full Markov)
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#' @param CPLBIdx Free variables correponding to matrix A_2
#' @param a "hyperrectangle" extension: bounds on hyperrectangle, "knapsack": objective function
#' @param ldelta Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param x.start Starting points for each chain. Needs to be fed in as a matrix with number of rows = number of chains
#' @param n.sample Number of samples
#' @param n.burnin Number of burnin steps
#' @param n.chains Number of chains
#' @param thinning Thinning parameter
#'
hyperrectangleSampler <- function(A,
                                  y,
                                  B = NULL,
                                  Model = "uniform",
                                  lambda = NULL,
                                  CPLBIdx = NULL,
                                  a = NULL,
                                  ldelta = 0,
                                  w = 0,
                                  x.start = NULL,
                                  n.sample = 1e+05,
                                  n.burnin = 1e+04,
                                  n.chains = 4,
                                  thinning = 1){
  r = nrow(A)
  c = ncol(A)

  if(round(n.sample/thinning) != n.sample/thinning){
    stop("thinning paramter must divide n.sample parameter!")
  }

  if(is.null(CPLBIdx)){
    CPLBIdx = (r+1):c
  } else if (length(CPLBIdx) != c-r){
    stop("CPLB index has to have length c-r")
  }

  if(is.null(B)){
    B = computeCPLB(A, CPLBIdx)
  }

  if(is.null(a)){
    a = numeric(c-r)

    for(i in 1:(c-r)){
      a[i] = lpSolve::lp(direction = "max",
                         objective.in = diag(c)[CPLBIdx[i],],
                         const.mat = A,
                         const.rhs = y,
                         const.dir = "=",
                         all.int = TRUE)$solution[CPLBIdx[i]]
    }
  } else if(length(a) == 1){
    a = rep(a, r)
  } else if(length(a) != c-r){
    stop("a has to have dimension c-r!")
  }

  cat("Initializing chains... \n")
  x = matrix(NA, ncol = c+3, nrow = n.chains*n.sample/thinning)
  x[, c+3] = 1:(n.chains*n.sample/thinning)
  x[, c+2] = rep(1:(n.sample/thinning), n.chains)
  x[, c+1] = rep(1:n.chains, each = n.sample/thinning)

  if(is.null(x.start)){
    x.start = matrix(NA, nrow = n.chains, ncol = c)
    for(i in 1:n.chains){
      x.start[i,] = lpSolve::lp(direction = "min",
                                objective.in = sample(0:1, c, replace = TRUE),
                                const.mat = A,
                                const.rhs = y,
                                const.dir = "=",
                                all.int = TRUE)$solution
    }
  }

  for(iii in 1:n.chains){
    x.current = x.start[iii,]

    moveIndices = sample(1:ncol(B), n.burnin, replace = TRUE)
    for(iiii in 1:n.burnin){
      moveIdx = moveIndices[iiii]
      x.proposal = proposePoint(x.current, moveIdx, B, extension = "hyperrectangle", CPLBIdx = CPLBIdx, a = a)
      alpha = acceptanceExt(x.current, x.proposal, ldelta, w, Model, lambda)

      if(stats::runif(1) < exp(alpha)){
        x.current = x.proposal
      }
    }

    moveIndices = sample(1:ncol(B), n.sample, replace = TRUE)
    for(iiiii in 1:n.sample){
      moveIdx = moveIndices[iiiii]
      x.proposal = proposePoint(x.current, moveIdx, B, extension = "hyperrectangle", CPLBIdx = CPLBIdx, a = a)
      alpha = acceptanceExt(x.current, x.proposal, ldelta, w, Model, lambda)

      if(stats::runif(1) < exp(alpha)){
        x.current = x.proposal
      }

      if(iiiii%%thinning == 0){
        x[(iii-1)*(n.sample/thinning) + iiiii/thinning, 1:c] = x.current
      }
    }
    cat("Chain ", iii, " completed.\n")
  }
  return(x)
}
