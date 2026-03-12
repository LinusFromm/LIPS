#' Extension Samplers
#' @param A Configuration Matrix
#' @param y (Corrupted/aggregated) Observation Vector
#' @param B Set of moves used (p-Markov, CPLB or Full Markov)
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#' @param extension Type of extension ("p", "hyperrectangle", "knapsack")
#' @param p Relaxation limit
#' @param CPLBIdx Free variables correponding to matrix A_2
#' @param a "hyperrectangle" extension: bounds on hyperrectangle, "knapsack": objective function
#' @param b Upper bound of objective function for knapsack extension
#' @param ldelta Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param x.start Starting points for each chain. Needs to be fed in as a matrix with number of rows = number of chains
#' @param n.sample Number of samples
#' @param n.burnin Number of burnin steps
#' @param n.chains Number of chains
#' @param thinning Thinning parameter
#'
#' @export
extensionSampler <- function(A,
                            y,
                            B = NULL,
                            Model = "uniform",
                            lambda = NULL,
                            extension = "p",
                            p = 0,
                            CPLBIdx = NULL,
                            a = NULL,
                            b = 0,
                            ldelta = 0,
                            w = 0,
                            x.start = NULL,
                            n.sample = 1e+05,
                            n.burnin = 1e+04,
                            n.chains = 4,
                            thinning = 1){

  if(extension == "p"){
    x = pSampler(A, y, B, Model, lambda, p, ldelta, w, x.start, n.sample, n.burnin, n.chains, thinning)
  } else if(extension == "hyperrectangle"){
    x = hyperrectangleSampler(A, y, B, Model, lambda, CPLBIdx, a, ldelta, w, x.start, n.sample, n.burnin, n.chains, thinning)
  } else if(extension == "knapsack"){
    x = knapsackSampler(A, y, B, Model, lambda, CPLBIdx, a, ldelta, w, x.start, n.sample, n.burnin, n.chains, thinning)
  }
  return(x)
}
