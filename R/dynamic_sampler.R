#' Dynamic lattice basis sampler
#' @param A Configuration matrix
#' @param y Observation vector
#' @param dist Distribution of the model ("Pois" - Poisson, "Unif" - Uniform), (Default = "Pois")
#' @param rate Poisson rate (if needed)
#' @param Method Sampling method ("MH" - Metropolis-Hasting, "Gibbs" - Gibbs sampling), (Default = "MH")
#' @param Proposal Move length proposal distribution ("Unif" - Uniform, "NonUnif" - Non-uniform), (Default = "Unif")
#' @param Reorder Should the columns of A be reordered? (Default = TRUE)
#' @param x_order_init If Reorder=FALSE, x_order can be used to reorder columns of A to match ordering of entries of x. Defaults to NULL when no such reordering is performed.
#' @param tune_par Tuning parameter for the dynamic Markov basis that controls variation of the fitness values for different CPLB (Default = 0.5)
#' @param B Set of moves used to move in the fibre. If NULL use dynamic sampler if not NULL use regular sampler. (Default = NULL)
#' @param x_init Initial state of the MCMC (not necessary, if nor provided algortihm will use lpSolve to find an initial solution)
#' @param nsample Number of samples per chain (Default = 100,000)
#' @param nburnin Number of burnin samples in each chain (Default = 10,000)
#' @param nchains Number of chains (Default = 4)
#' @param thinning Record only every kth sample where k is the value for thinning (Default = 1)
#' @param combine Should extra moves be included combining lattice basis vectors? Defaults to FALSE, but should usually be set to TRUE if A is not unimodular.
#' @param verbose Controls level of detail in recording lattice bases used. (Namely, if TRUE it records the x_order at every iteration)
#'
#' @return A list with components X (a matrix, each row corresponding to samplers for an entry of x) and x_order (a vector describing dynamic selection of lattice bases, if verbose=1).
#' @export
dynamic_sampler <- function(A = NULL,
                            y = NULL,
                            dist = "Pois",
                            rate = NULL,
                            Method = "MH",
                            Proposal = "Unif",
                            Reorder = TRUE,
                            x_order_init = NULL,
                            tune_par = 0.5,
                            B = NULL,
                            x_init = NULL,
                            nsample = 1e+05,
                            nburnin = 1e+04,
                            nchains = 4,
                            thinning = 1,
                            combine = FALSE,
                            verbose = FALSE){

  # If dist = "Unif" ensure we only use "Gibbs" sampling and never "MH"
  if(dist == "Unif" && Method == "MH"){
    warning("Only method for the uniform model is Gibbs!")
    Method = "Gibbs"
  }

  # Why does the proposal have to be uniform now?
  if(combine && Proposal == "NonUnif"){
    warning("Proposal must be uniform if combine == TRUE")
    Proposal = "Unif"
  }

  y = as.numeric(y)

  ## Delete all zero columns
  zero_col_idx = apply(A, 2, function(col) all(col == 0))
  non_zero_col_idx = !zero_col_idx
  zero_cols = sum(zero_col_idx)

  if(zero_cols > 0){
    # Do I need to store c_full and rate_full?
    c_full = ncol(A)
    rate_full = rate

    A = A[,non_zero_col_idx]
    rate = rate[non_zero_col_idx]
    x_init = x_init[non_zero_col_idx]
    x_order = x_order[non_zero_col_idx]
  }

  c = ncol(A)
  r = nrow(A)
  tol = 1e-19

  # If a Markov basis is provided we do not need to reorder x and no dynamic sampling
  if(!is.null(B)){
    cat("\nDynamic sampling is deactivated by providing a Markov basis...")
    Reorder = FALSE
    dA1 = 1
    tune_par = -1
  }
  # Matrix can be automatically ordered so that the first n columns are linearly independent (Reorder == TRUE)
  # or we use a manually provided x_order to input an order that has A with first n columns linearly independent (Reorder == FALSE && x_order != NULL)
  if(Reorder){
    rate_star = rate
    lambda_order = order(rate_star, decreasing = FALSE)
    A = A[, lambda_order]
    x_order_init = lambda_order
    for(i in 3:r){
      while(qr(A[,1:i])$rank < i){
        A = A[,c(1:(i-1),(i+1):c,i)]
        x_order_init = x_order_init[c(1:(i-1),(i+1):c,i)]
      }
    }
  } else {
    if(is.null(x_order_init)) x_order_init = 1:c
    A = A[, x_order_init]
  }

  rate_init = rate[x_order_init]

  # If no B provided we produce the first CPLB
  if(is.null(B)){
    A1 = A[,1:r]
    A2 = A[,-(1:r)]
    dA1 = det(A1)
    C_init = solve(A1)%*%A2
    B_init = rbind(-solve(A1)%*%A2, diag(c-r))

    # If B is non-integer then the moves must be rescaled. We do this later during the actual sampling
  }

  cat("\nInitializing chains...")
  m = c-r+1*combine
  X = matrix(0, ncol = c + 3, nrow = nchains*nsample)

  if(verbose){
    X_ORDER = matrix(0, ncol = c+3, nrow = nchains*nsample)
  }

  for(chain in 1:nchains){
    x_order = x_order_init
    rate = rate_init

    C = C_init
    B = B_init
    if(is.null(x_init)){
      x_init <- lpSolve::lp("max",
                            objective.in=sample(0:1, c, replace = TRUE),
                            const.mat=A,
                            const.dir=rep("=",nrow(A)),
                            const.rhs=y,
                            all.int=TRUE)$solution
    }
    x = x_init

    for (burin in nburnin) {
      j = sample(1:m, 1)

      # If sampled j outside of range of 1:ncol(B) we generate a move from moves in B. (This makes this algorithm irreducible when CPLBs are not a Markov basis for the fibre)
      if(j <= ncol(B)) {
        z = B[,j]
      } else {
        delta = 0
        while(sum(delta) == 0) delta = stats::rpois(ncol(B), lambda = 0.5)
        delta = delta*sample(c(-1,1), size = ncol(B), replace = TRUE)
        z = B%*%delta
      }

      if(abs(dA1) != 1){
        # A is not unimodular we need to rescale the basis moves
        if(is_wholenumber(z) == FALSE) z = round(z*abs(dA1))/numbers::mGCD(round(abs(z*dA1)))
      }
      max_move = floor(min((x/abs(z))[which(z<0)]))
      min_move = -floor(min((x/abs(z))[which(z>0)]))
      x_min = x+min_move*z
      x_max = x+max_move*z
      idx = 1:(max_move-min_move+1)

      #Slight speed-up by only updating the entries that are non-zero in z
      update_idx = which(z!=0)
      if(max(idx > 1)){
        if(Method == "Gibbs"){
          if(dist != "Unif") {
             x_matrix = round(t(mapply(seq,
                                       from=x_min[update_idx],
                                       by = z[update_idx],
                                       length.out = max_move-min_move+1)))
          }
          if(dist == "Pois"){
             log_probs = colSums(stats::dpois(x = x_matrix,
                                              lambda = rate[update_idx],
                                              log = TRUE))

             probs = exp(log_probs - max(log_probs))
             x = x_min + sample(x = idx-1, size = 1, prob = probs)*z
          }
          if(dist == "Unif"){
            x = x_min + sample(idx-1, size = 1)*z
          }
        } else if(Method == "MH"){
          if(Proposal == "Unif"){
            move_length = sample(min_move:max_move, 1)
          } else if(Method == "NonUnif"){
            if(dist == "Pois"){
              aa = x[r+j] + z[r+j]*min_move
              bb = x[r+j] + z[r+j]*max_move
              move_length = (extraDistr::rtpois(1, lambda=rate[r+j], a = aa-0.5, b = bb)-x[r+j])/z[r+j]
            }
          }
          x_cand = x+z*move_length
          if(dist == "Pois"){
            L = sum(stats::dpois(x[update_idx],
                    lambda = rate[update_idx],
                    log = TRUE))
            L_cand = sum(stats::dpois(x_cand[update_idx],
                    lambda = rate[update_idx],
                    log = TRUE))
          }
          if(Proposal == "Unif"){
            acc_prob = exp(L_cand-L)
          } else if(Proposal == "NonUnif"){
            if(dist == "Pois"){
              q_can = stats::dpois(x_cand[r+j],
                                   lambda = rate[r+j],
                                   log = TRUE)
              q_cur = stats::dpois(x[r+j],
                                   lambda = rate[r+j],
                                   log = TRUE)
            }
            acc_prob = exp(L_cand - L + q_cur - q_can)
          }
          if(is.na(acc_prob)) acc_prob = 0
          if(stats::runif(1) < acc_prob) x = x_cand
        }
      }

      if(tune_par > 0){
        rate_star = stats::rnorm(c, mean = rate, sd = tune_par*rate)
        ii = sample(1:r,1)
        swap_idx = abs(C[ii,])>tol

        if(any(swap_idx)){
          jj = sample((1:(c-r))[swap_idx],1)
          if(rate_star[ii] <= rate_star[jj+r]/abs(C[ii,jj])){
            ei = rep(0,r)
            ej = rep(0,c-r)
            ei[ii] = 1
            ej[jj] = 1
            dA1 = round(C[ii,jj]*dA1)
            C = C-outer(C[,jj]-ei, C[ii,]+ej)/C[ii,jj]
            B = rbind(-C, diag(c-r))
            x_order[c(ii,jj+r)] = x_order[c(jj+r,ii)]
            x[c(ii,jj+r)] = x[c(jj+r,ii)]
            rate[c(ii, jj+r)] = rate[c(jj+r,ii)]
          }
        }
      }
      x = round(x,10)
    }

    for (iter in 1:nsample) {
      j = sample(1:m, 1)

      # If sampled j outside of range of 1:ncol(B) we generate a move from moves in B. (This makes this algorithm irreducible when CPLBs are not a Markov basis for the fibre)
      if(j <= ncol(B)) {
        z = B[,j]
      } else {
        delta = 0
        while(sum(delta) == 0) delta = stats::rpois(ncol(B), lambda = 0.5)
        delta = delta*sample(c(-1,1), size = ncol(B), replace = TRUE)
        z = B%*%delta
      }

      if(abs(dA1) != 1){
        # A is not unimodular we need to rescale the basis moves
        if(is_wholenumber(z) == FALSE) z = round(z*abs(dA1))/numbers::mGCD(round(abs(z*dA1)))
      }
      max_move <- floor(min((x/abs(z))[which(z<0)]))
      min_move <- -floor(min((x/abs(z))[which(z>0)]))
      x_min <- x+min_move*z
      x_max <- x+max_move*z
      idx <- 1:(max_move-min_move+1)

      #Slight speed-up by only updating the entries that are non-zero in z
      update_idx = which(z!=0)
      if(max(idx > 1)){
        if(Method == "Gibbs"){
          if(dist != "Unif") {
            x_matrix = round(t(mapply(seq,
                                      from=x_min[update_idx],
                                      by = z[update_idx],
                                      length.out = max_move-min_move+1)))
          }
          if(dist == "Pois"){
            log_probs = colSums(stats::dpois(x = x_matrix,
                                             lambda = rate[update_idx],
                                             log = TRUE))

            probs = exp(log_probs - max(log_probs))
            x = x_min + sample(x = idx-1, size = 1, prob = probs)*z
          }
          if(dist == "Unif"){
            x = x_min + sample(idx-1, size = 1)*z
          }
        } else if(Method == "MH"){
          if(Proposal == "Unif"){
            move_length = sample(min_move:max_move, 1)
          } else if(Method == "NonUnif"){
            if(dist == "Pois"){
              aa = x[r+j] + z[r+j]*min_move
              bb = x[r+j] + z[r+j]*max_move
              move_length = (extraDistr::rtpois(1, lambda=rate[r+j], a = aa-0.5, b = bb)-x[r+j])/z[r+j]
            }
          }
          x_cand = x+z*move_length
          if(dist == "Pois"){
            L = sum(stats::dpois(x[update_idx],
                    lambda = rate[update_idx],
                    log = TRUE))
            L_cand = sum(stats::dpois(x_cand[update_idx],
                         lambda = rate[update_idx],
                         log = TRUE))
          }
          if(Proposal == "Unif"){
            acc_prob = exp(L_cand-L)
          } else if(Proposal == "NonUnif"){
            if(dist == "Pois"){
              q_can = stats::dpois(x_cand[r+j],
                                   lambda = rate[r+j],
                                   log = TRUE)
              q_cur = stats::dpois(x[r+j],
                                   lambda = rate[r+j],
                                   log = TRUE)
            }
            acc_prob = exp(L_cand - L + q_cur - q_can)
          }
          if(is.na(acc_prob)) acc_prob = 0
          if(stats::runif(1) < acc_prob) x = x_cand
        }
      }

      if(tune_par > 0){
        rate_star = stats::rnorm(c, mean = rate, sd = tune_par*rate)
        ii = sample(1:r,1)
        swap_idx = abs(C[ii,])>tol

        if(any(swap_idx)){
          jj = sample((1:(c-r))[swap_idx],1)
          if(rate_star[ii] <= rate_star[jj+r]/abs(C[ii,jj])){
            ei = rep(0,r)
            ej = rep(0,c-r)
            ei[ii] = 1
            ej[jj] = 1
            dA1 = round(C[ii,jj]*dA1)
            C = C-outer(C[,jj]-ei, C[ii,]+ej)/C[ii,jj]
            B = rbind(-C, diag(c-r))
            x_order[c(ii,jj+r)] = x_order[c(jj+r,ii)]
            x[c(ii,jj+r)] = x[c(jj+r,ii)]
            rate[c(ii, jj+r)] = rate[c(jj+r,ii)]
          }
        }
      }

      x = round(x,10)
      X[(chain-1)*nsample+iter, x_order] <- x
      X[(chain-1)*nsample+iter, c + (1:3)] <- c(chain, iter, (chain-1)*nsample+iter)

      if (verbose==1) {
        X_ORDER[(chain-1)*nsample+iter, 1:c] <- x_order
        X_ORDER[(chain-1)*nsample+iter, c+(1:3)] <- c(chain, iter, (chain-1)*nsample+iter)
      }
    }

    x_init = NULL
    cat("\nChain", chain, "completed!")
  }

  if(verbose == 1) x_order = X_ORDER
  if(zero_cols > 0){
    XX = X
    X = matrix(NA, ncol = c_full+3, nrow = nrow(XX))
    X[, c(non_zero_col_idx, (ncol(X)-2:ncol(X)))] = XX
    for (i in 1:zero_cols){
      X[,which(zero_col_idx)[i]] <- stats::rpois(nrow(XX), rate_full[which(zero_col_idx)[i]])
    }
  }

  return(list(x = X, x_order = x_order))
}
