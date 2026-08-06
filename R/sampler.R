#' Linear inverse problem Samplers
#' @param A Configuration Matrix
#' @param y (Corrupted/aggregated) Observation Vector
#' @param B Set of moves used (p-Markov, CPLB or Full Markov)
#' @param Model Type of model ("uniform" or "poisson")
#' @param lambda Rate parameters for "poisson" Model
#' @param extension Type of extension ("p", "hyperrectangle", "knapsack")
#' @param p Relaxation limit
#' @param A2Idx Free variables correponding to matrix A_2
#' @param a "hyperrectangle" extension: bounds on hyperrectangle, "knapsack": objective function
#' @param b Upper bound of objective function for knapsack extension
#' @param loggamma Combined normalizing constant and delta
#' @param w Weight put on "outsideness"
#' @param xStart Starting points for each chain. Needs to be fed in as a matrix with number of rows = number of chains
#' @param nSample Number of samples
#' @param nBurnin Number of burnin steps
#' @param nChains Number of chains
#' @param thinning Thinning parameter
#' @param includeBurnin Should the output include the burn-in samples?
#'
#' @export
sampler <- function(A,
                    y,
                    B = NULL,
                    Model = "uniform",
                    lambda = NULL,
                    extension = "none",
                    p = 0,
                    A2Idx = NULL,
                    a = NULL,
                    b = 0,
                    loggamma = 0,
                    w = 0,
                    xStart = NULL,
                    nSample = 1e+05,
                    nBurnin = 1e+04,
                    nChains = 4,
                    thinning = 1,
                    includeBurnin = FALSE){

  future::plan(future::multisession, workers = min(nChains, future::availableCores()))

  if(extension == "none"){
    x = future.apply::future_lapply(1:nChains, function(chainID) {
      xSampler(
        A = A,
        y = y,
        B = B,
        Model = Model,
        lambda = lambda,
        xStart = xStart,
        nSample = nSample,
        nBurnin = nBurnin,
        chainID = chainID,
        thinning = thinning,
        includeBurnin = includeBurnin
      )
    }, future.seed = TRUE, future.packages = "LIPS")
  }

  if(extension == "p"){
    x = future.apply::future_lapply(1:nChains, function(chainID) {
      pSampler(
        A = A,
        y = y,
        B = B,
        Model = Model,
        lambda = lambda,
        p = p,
        loggamma = loggamma,
        w = w,
        xStart = xStart,
        nSample = nSample,
        nBurnin = nBurnin,
        chainID = chainID,
        thinning = thinning,
        includeBurnin = includeBurnin
      )
    }, future.seed = TRUE, future.packages = "LIPS")
  } else if(extension == "hyperrectangle"){
    x = future.apply::future_lapply(1:nChains, function(chainID) {
      hyperrectangleSampler(A = A,
                            y = y,
                            B = B,
                            Model = Model,
                            lambda = lambda,
                            A2Idx = A2Idx,
                            a = a,
                            loggamma = loggamma,
                            w = w,
                            xStart = xStart,
                            nSample = nSample,
                            nBurnin = nBurnin,
                            chainID = chainID,
                            thinning = thinning)
    }, future.seed = TRUE, future.packages = "LIPS")
  } else if(extension == "knapsack"){
    x = future.apply::future_lapply(1:nChains, function(chainID) {
      knapsackSampler(A = A,
                      y = y,
                      B = B,
                      Model = Model,
                      lambda = lambda,
                      A2Idx = A2Idx,
                      a = a,
                      b = b,
                      loggamma = loggamma,
                      w = w,
                      xStart = xStart,
                      nSample = nSample,
                      nBurnin = nBurnin,
                      chainID = chainID,
                      thinning = thinning)
    }, future.seed = TRUE, future.packages = "LIPS")
  }
  return(do.call(rbind, x))
}
