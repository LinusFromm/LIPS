#' Compute a specific fibre
#' @param A Configuration matrix
#' @param y Observation vector
#' @param epsilon Tolerance for off-integerness
#'
#' @return Returns a Matrix containing all points of the fibre as columns.
#' @export
generate_fibre <- function(A,
                           y,
                           epsilon = 1e-04){
  x_max = c()
  x_min = c()

  for(i in 1:ncol(A)){
    x_max[i] = lpSolve::lp(direction = "max",
                           objective.in = diag(ncol(A))[i,],
                           const.mat = A,
                           const.dir = "=",
                           const.rhs = y,
                           all.int = TRUE)$solution[i]
    x_min[i] = lpSolve::lp(direction = "min",
                           objective.in = diag(ncol(A))[i,],
                           const.mat = A,
                           const.dir = "=",
                           const.rhs = y,
                           all.int = TRUE)$solution[i]
  }

  idx = c(1)
  i = 2
  while(qr(A[,idx])$rank < qr(A)$rank){
    temp = A[,c(idx,i)]
    if(qr(temp)$rank > qr(A[,idx])$rank){
      idx = c(idx, i)
    }
    i = i+1
  }

  hyperrectangle_ext_coord = setdiff(1:ncol(A), idx)
  coords = apply(rbind(x_max, x_min)[,hyperrectangle_ext_coord], 2, function(col) return((col[2]):(col[1])), simplify = FALSE)

  grid = expand.grid(coords)
  A1 = A[,idx]
  A2 = A[,hyperrectangle_ext_coord]

  ext = c(solve(A1)%*%y)-solve(A1)%*%A2%*%t(grid)
  fibre_idx = which(apply(ext, 2, function(col) all(round(col) >= 0) && all(abs(col - round(col)) <= epsilon)))

  fibre = matrix(0, ncol = length(fibre_idx), nrow = ncol(A))
  fibre[idx,] = round(ext[,fibre_idx])
  fibre[hyperrectangle_ext_coord,] = round(t(grid)[,fibre_idx])

  return(fibre)
}
