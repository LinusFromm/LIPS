#' Plot traceplots for the variables
#' @param x The output from the MCMC should be put in.
#' @param variables Select a subset of variables whose traceplot should be plotted
#'
#' @return Plots a traceplot for each variable in the subset defined in variables
#' @export
traceplot <- function(x,
                      variables){
  n = ncol(x)
  m = nrow(x)

  n.chains = max(x[,n-2])
  n.samples = max(x[,n-1])

  for(i in 1:length(variables)){
    plot(x[1:n.samples,variables[i]],
         type = "l",
         col = scales::alpha(1, 1/n.chains),
         ylab = paste0("x_",variables[i]),
         xlab = "Iteration",
         ylim = c(min(x[,variables[i]]), max(x[,variables[i]])))

    for(j in 1:(n.chains-1)){
      graphics::lines(x[(j*n.samples+1):((j+1)*n.samples), variables[i]], col = scales::alpha(j+1, 1/n.chains))
    }
  }
}
