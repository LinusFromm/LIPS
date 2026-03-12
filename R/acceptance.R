acceptance <- function(x, x.star, Model, lambda = NULL){
  x <- as.numeric(x)
  x.star <- as.numeric(x.star)
  lambda <- as.numeric(lambda)

  if(Model == "poisson"){
    alpha = min(0, sum(stats::dpois(x.star, lambda = lambda, log = TRUE)) - sum(stats::dpois(x, lambda = lambda, log = TRUE)))
  }
  return(alpha)
}
