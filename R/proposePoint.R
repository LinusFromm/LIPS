#' Proposing algorithm
#' @param x.current Current state of MCMC
#' @param moveIdx Index of the move in the collection of moves B
#' @param B Collection of moves (p-MB, CPLB or FMB)
#' @param extension Type of extension
#' @param p Relaxation limit
#' @param A2Idx Free variables correponding to matrix A_2
#' @param a "hyperrectangle" extension: bounds on hyperrectangle, "knapsack": objective function
#' @param b Upper bound of objective function for knapsack extension
#'
proposePoint <- function(x.current, moveIdx, B, extension = "p", p = 0, A2Idx = NULL, a = NULL, b = NULL){
  z = B[, moveIdx]

  if(extension == "p"){

    cmax = floor(min(((x.current+p)/abs(z))[which(z<0)]))
    cmin = -floor(min(((x.current+p)/abs(z))[which(z>0)]))

    xmin = x.current + cmin*z
    c = sample(0:(cmax-cmin), 1, replace = TRUE)

  } else if(extension == "hyperrectangle"){

    xmin = x.current - x.current[A2Idx[moveIdx]]*z
    c = sample(0:a[moveIdx], 1, replace = TRUE)

  } else if(extension == "knapsack"){

    xmin = x.current - x.current[A2Idx[moveIdx]]*z
    cmax = floor((b-sum(xmin[A2Idx]*a))/a[moveIdx])
    c = sample(0:cmax, 1, replace = TRUE)

  }
  return(xmin + c*z)
}
