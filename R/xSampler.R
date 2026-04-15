#' x-sampler
#' @param A Configuration matrix
#' @param y Vector of observations
#' @param B Full Markov bases to be used
#' @param Model Unconditional distribution on x
#' @param lambda Poisson rates for each x_i
#' @param proposal Which proposal for moves? "Random" or "NonRandom", first will use the uniform distribution the latter a systematic approach.
#' @param nSample Number of samples that are produced
#' @param nChains Number of chains the sample is produced from
#' @param nBurnin Number of burnin steps
#' @param thinning Thinning parameter (store only every ith state)
#'
#' @return Returns a matrix of samples
#' @export
xSampler <- function(A,
                      y,
                      B,
                      Model = "Unif",
                      lambda = NULL,
                      proposal = "Random",
                      nSample = 1e+05,
                      nBurnin = 1e+04,
                      nChains = 4,
                      thinning = 1){
  c = ncol(A)
  r = nrow(A)
  d = c-r
  m = ncol(B)

  cat("Initializing chains... \n")
  x = matrix(0, nrow = nChains*(nSample/thinning), ncol = c+3)

  for(ii in 1:nChains){
    xCurrent = lpSolve::lp("max",
                            objective.in = sample(0:1, c, replace = TRUE),
                            const.mat = A,
                            const.dir = "==",
                            const.rhs = y,
                            all.int = TRUE)$solution

    move_indices = c()
    if(proposal == "Random"){
      move_indices = sample.int(m, nBurnin, replace = TRUE)
    } else if (proposal == "NonRandom"){
      move_indices = c(rep(1:m, floor(nBurnin/m)), 1:(nBurnin%%m))
    }
    for(iii in 1:nBurnin){
      z.idx = move_indices[iii]
      z = B[,z.idx]

      amax = floor(min((xCurrent/abs(z))[which(z<0)]))
      amin = -floor(min((xCurrent/abs(z))[which(z>0)]))

      # Try uniform first
      xmin = xCurrent + amin*z
      a = sample(0:(amax-amin), 1, replace = TRUE)

      # Here I need to also use the MH-acceptance ratio when implementing other distributions
      if(Model == "Unif"){
        xCurrent = xmin + a*z
      } else if(Model == "Pois") {
        alpha = sum(stats::dpois(xmin + a*z, lambda = lambda, log = TRUE) - stats::dpois(xCurrent, lambda = lambda, log = TRUE))
        u = stats::runif(1, 0, 1)

        if(u < exp(alpha)){
          xCurrent = xmin + a*z
        }
      }
    }

    move_indices = c()
    if(proposal == "Random"){
      move_indices = sample.int(m, nSample, replace = TRUE)
    } else if (proposal == "NonRandom"){
      move_indices = c(rep(1:m, floor(nSample)/m), 1:(nSample%%m))
    }
    for(iiii in 1:nSample){
      z.idx = move_indices[iiii]
      z = B[,z.idx]

      amax = floor(min((xCurrent/abs(z))[which(z<0)]))
      amin = -floor(min((xCurrent/abs(z))[which(z>0)]))

      # Try uniform first
      xmin = xCurrent + amin*z
      a = sample(0:(amax-amin), 1, replace = TRUE)

      ## Here I need to also use the MH-acceptance ratio when implementing other distributions
      if(Model == "Unif"){
        xCurrent = xmin + a*z
      } else if(Model == "Pois") {
        alpha = sum(stats::dpois(xmin + a*z, lambda = lambda, log = TRUE) - stats::dpois(xCurrent, lambda = lambda, log = TRUE))
        u = stats::runif(1, 0, 1)

        if(u < exp(alpha)){
          xCurrent = xmin + a*z
        }
      }

      if(iiii %% thinning == 0){
        x[(ii-1)*(nSample/thinning) + iiii/thinning, 1:c] = xCurrent
        x[(ii-1)*(nSample/thinning) + iiii/thinning, (c+1):(c+3)] = c(ii, iiii/thinning, (ii-1)*nSample/thinning + iiii/thinning)
      }
    }
    cat("Chain ", ii, " completed.\n")
  }

  return(x)
}
