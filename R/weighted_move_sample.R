#' Weighted moves sampler
#' @param A Configuration matrix
#' @param y Vector of observations
#' @param B Full Markov bases to be used
#' @param dist Unconditional distribution on x
#' @param rate Poisson rates for each x_i
#' @param n.sample Number of samples that are produced
#' @param n.chains Number of chains the sample is produced from
#' @param n.burnin Number of burnin steps
#' @param thinning Thinning parameter (store only every ith state)
#' @param a Tuning parameter determining how much heavier weights count towards sampling of moves
#' @param weights Weights for the individual moves. If NULL uses uniform weights. (Default - NULL)
#'
#' @return Returns a matrix of samples
#' @export
weighted_move_sampler <- function(A,
                      y,
                      B,
                      dist = "Unif",
                      rate = NULL,
                      n.sample = 1e+05,
                      n.burnin = 1e+04,
                      n.chains = 4,
                      thinning = 1,
                      a = 1,
                      weights = NULL){
  c = ncol(A)
  r = nrow(A)
  d = c-r
  m = ncol(B)
  if(is.null(weights)){
    weights = rep(1, m)
  }

  cat("Initializing chains... \n")
  x = matrix(0, nrow = n.chains*(n.sample/thinning), ncol = c+3)



  for(ii in 1:n.chains){
    x.current = lpSolve::lp("max",
                            objective.in = sample(0:1, c, replace = TRUE),
                            const.mat = A,
                            const.dir = "==",
                            const.rhs = y,
                            all.int = TRUE)$solution

    move_indices = sample.int(m, n.burnin, replace = TRUE, prob = weights^a)
    for(iii in 1:n.burnin){
      z.idx = move_indices[iii]
      z = B[,z.idx]

      amax = floor(min((x.current/abs(z))[which(z<0)]))
      amin = -floor(min((x.current/abs(z))[which(z>0)]))

      # Try uniform first
      xmin = x.current + amin*z
      a = sample(0:(amax-amin), 1, replace = TRUE)

      # Here I need to also use the MH-acceptance ratio when implementing other distributions
      if(dist == "Unif"){
        x.current = xmin + a*z
      } else if(dist == "Pois") {
        alpha = sum(stats::dpois(xmin + a*z, lambda = rate, log = TRUE) - stats::dpois(x.current, lambda = rate, log = TRUE))
        u = stats::runif(1, 0, 1)

        if(u < exp(alpha)){
          x.current = xmin + a*z
        }
      }
    }

    move_indices = sample.int(m, n.sample, replace = TRUE, prob = weights^a)
    for(iiii in 1:n.sample){
      z.idx = move_indices[iiii]
      z = B[,z.idx]

      amax = floor(min((x.current/abs(z))[which(z<0)]))
      amin = -floor(min((x.current/abs(z))[which(z>0)]))

      # Try uniform first
      xmin = x.current + amin*z
      a = sample(0:(amax-amin), 1, replace = TRUE)

      ## Here I need to also use the MH-acceptance ratio when implementing other distributions
      if(dist == "Unif"){
        x.current = xmin + a*z
      } else if(dist == "Pois") {
        alpha = sum(stats::dpois(xmin + a*z, lambda = rate, log = TRUE) - stats::dpois(x.current, lambda = rate, log = TRUE))
        u = stats::runif(1, 0, 1)

        if(u < exp(alpha)){
          x.current = xmin + a*z
        }
      }

      if(iiii %% thinning == 0){
        x[(ii-1)*(n.sample/thinning) + iiii/thinning, 1:c] = x.current
        x[(ii-1)*(n.sample/thinning) + iiii/thinning, (c+1):(c+3)] = c(ii, iiii/thinning, (ii-1)*n.sample/thinning + iiii/thinning)
      }
    }
    cat("Chain ", ii, " completed.\n")
  }

  return(x)
}
