#' Converting from matrix to draws_df
#' @param x Matrix containing the sample including chain, iteration and draw numbers.
#' @param var.names Fill names for the variables if NULL names will be $xi$ for the ith variable.
#' @param matrix Is the sample a matrix?
#'
#' @return Returns the sample object of any of the samplers as a draws_df file
#' @export
to_draws_df <- function(x,
                        var.names = NULL,
                        matrix = TRUE){

  if(matrix){
    x_df = as.data.frame(x)
  } else {
    x_df = as.data.frame(do.call(rbind, x))
  }

  if(is.null(var.names)){
    colnames(x_df) <- c(paste0("x",1:(ncol(x_df)-3)), ".chain", ".iteration", ".draw")
  } else {
    colnames(x_df) <- c(var.names, ".chain", ".iteration", ".draw")
  }
  return(posterior::as_draws_df(x_df))
}
