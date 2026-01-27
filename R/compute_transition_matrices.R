#' Computes transition matrices for all moves
#' @param fibre A matrix that contains all every point as column
#' @param B A matrix that contains everey move of a basis as columns
#' @param epsilon Tolerance to off-integerness
#'
#' @return Returns a list of matrices. Each is the transition matrix for a particular move of B. These can be used to find the overall transition matrix for the whole setup.
#' @export
compute_transition_matrices <- function(fibre,
                                        B,
                                        epsilon = 1e-04){
  H.z = list()
  for(k in 1:ncol(B)){
    b = B[,k]
    H = matrix(FALSE, ncol = ncol(fibre), nrow = ncol(fibre))
    for(i in 1:ncol(fibre)){
      for(j in i:ncol(fibre)){
        connection = fibre[,i]-fibre[,j]
        non_zero_idx = which(connection!=0)

        H[i,j] = (length(unique(connection[non_zero_idx]/b[non_zero_idx]))==1) && all(abs(b[setdiff(1:nrow(fibre), non_zero_idx)]) <= epsilon)
      }
    }
    H = (H+t(H)+diag(ncol(fibre)))/rowSums(H+t(H)+diag(ncol(fibre)))
    H.z[[k]] = H
  }

  return(H.z)
}
