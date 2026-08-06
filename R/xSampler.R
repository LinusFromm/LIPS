#' x-sampler
#' @param A Configuration matrix
#' @param y Vector of observations
#' @param B Full Markov bases to be used
#' @param Model Unconditional distribution on x
#' @param lambda Poisson rates for each x_i
#' @param proposal Which proposal for moves? "Random" or "NonRandom", first will use the uniform distribution the latter a systematic approach.
#' @param xStart Starting points for each chain. Needs to be fed in as a matrix with number of rows = number of chains
#' @param nSample Number of samples that are produced
#' @param chainID Chain ID passed down from sampler.R
#' @param nBurnin Number of burnin steps
#' @param thinning Thinning parameter (store only every ith state)
#' @param includeBurnin Should the output include the burn-in samples?
#'
#' @return Returns a matrix of samples
#' @export
xSampler <- function(A,
                     y,
                     B,
                     Model = "uniform",
                     lambda = NULL,
                     proposal = "Random",
                     xStart = NULL,
                     nSample = 1e+05,
                     nBurnin = 1e+04,
                     chainID = 1,
                     thinning = 1,
                     includeBurnin = FALSE){
  c = ncol(A)
  r = nrow(A)
  d = c-r
  m = ncol(B)

  if(includeBurnin){
    x = matrix(NA, ncol = c+3, nrow = nSample/thinning+nBurnin/thinning)
    x[, c+3] = (chainID-1)*((nSample + nBurnin)/thinning) + 1:((nSample + nBurnin)/thinning)
    x[, c+2] = 1:((nSample + nBurnin)/thinning)
    x[, c+1] = rep(chainID, (nSample + nBurnin)/thinning)
  } else {
    x = matrix(NA, ncol = c+3, nrow = nSample/thinning)
    x[, c+3] = (chainID-1)*(nSample/thinning) + 1:(nSample/thinning)
    x[, c+2] = 1:(nSample/thinning)
    x[, c+1] = rep(chainID, nSample/thinning)
  }

  if(is.null(xStart)){
    xStart = lpSolve::lp(direction = "min",
                         objective.in = sample(0:1, c, replace = TRUE),
                         const.mat = A,
                         const.rhs = y,
                         const.dir = "=",
                         all.int = TRUE)$solution
  }

  xCurrent = xStart

  if(nBurnin > 0){
    temperature = seq(0,1,length.out = nBurnin)

    move_indices = c()
    if(proposal == "Random"){
      move_indices = sample.int(m, nBurnin, replace = TRUE)
    } else if (proposal == "NonRandom"){
      move_indices = c(rep(1:m, floor(nBurnin/m)), 1:(nBurnin%%m))
    }

    if(includeBurnin){
      for(iiii in 1:nBurnin){
        z.idx = move_indices[iiii]
        z = B[,z.idx]

        amax = floor(min((xCurrent/abs(z))[which(z<0)]))
        amin = -floor(min((xCurrent/abs(z))[which(z>0)]))

        # Try uniform first
        xmin = xCurrent + amin*z
        a = sample(0:(amax-amin), 1, replace = TRUE)

        # Here I need to also use the MH-acceptance ratio when implementing other distributions
        if(Model == "uniform"){
          xCurrent = xmin + a*z
        } else if(Model == "poisson") {
          alpha = sum(stats::dpois(xmin + a*z, lambda = lambda, log = TRUE) - stats::dpois(xCurrent, lambda = lambda, log = TRUE))
          u = stats::runif(1, 0, 1)

          if(u < exp(alpha*temperature[iiii])){
            xCurrent = xmin + a*z
          }

          if(iiii%%thinning == 0){
            x[iiii/thinning, 1:c] <- xCurrent
          }
        }
      }

    } else {
      for(iiii in 1:nBurnin){
        z.idx = move_indices[iiii]
        z = B[,z.idx]

        amax = floor(min((xCurrent/abs(z))[which(z<0)]))
        amin = -floor(min((xCurrent/abs(z))[which(z>0)]))

        # Try uniform first
        xmin = xCurrent + amin*z
        a = sample(0:(amax-amin), 1, replace = TRUE)

        # Here I need to also use the MH-acceptance ratio when implementing other distributions
        if(Model == "uniform"){
          xCurrent = xmin + a*z
        } else if(Model == "poisson") {
          alpha = sum(stats::dpois(xmin + a*z, lambda = lambda, log = TRUE) - stats::dpois(xCurrent, lambda = lambda, log = TRUE))
          u = stats::runif(1, 0, 1)

          if(u < exp(alpha*temperature[iiii])){
            xCurrent = xmin + a*z
          }
        }
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
    if(Model == "uniform"){
      xCurrent = xmin + a*z
    } else if(Model == "poisson") {
      alpha = sum(stats::dpois(xmin + a*z, lambda = lambda, log = TRUE) - stats::dpois(xCurrent, lambda = lambda, log = TRUE))
      u = stats::runif(1, 0, 1)

      if(u < exp(alpha)){
        xCurrent = xmin + a*z
      }
    }

    if(iiii %% thinning == 0){
      if(includeBurnin){
        x[nBurnin/thinning + iiii/thinning, 1:c] = xCurrent
      } else {
        x[iiii/thinning, 1:c] = xCurrent
      }
    }
  }
  return(x)
}
