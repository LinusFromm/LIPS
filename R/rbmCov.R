#' Compute the multidimensional ESS for the sample
#' @param chains Sampled chains in list format (each item is a chain in matrix format without indices)
#' @param batchSize Variable used to evaluate the covarince matrix
#'
rbmCov <- function(chains, batchSize = NULL) {
  chains <- lapply(chains, as.matrix)

  m <- length(chains)
  n <- unique(vapply(chains, nrow, integer(1)))
  d <- unique(vapply(chains, ncol, integer(1)))

  if (length(n) != 1L) {
    stop("This implementation requires equal chain lengths.")
  }

  if (length(d) != 1L) {
    stop("All chains must have the same number of columns.")
  }

  if (is.null(batchSize)) {
    batchSize <- floor(sqrt(n))
  }

  b <- batchSize
  a <- floor(n / b)
  n_used <- a * b

  if (m * a <= d) {
    warning(
      "Too few total batches relative to the dimension; ",
      "the covariance estimate may be singular."
    )
  }

  # Discard at most b - 1 final observations from each chain
  chains_used <- lapply(
    chains,
    function(z) z[seq_len(n_used), , drop = FALSE]
  )

  global_mean <- colMeans(do.call(rbind, chains_used))

  batch_means <- do.call(
    rbind,
    lapply(chains_used, function(z) {
      batch_id <- rep(seq_len(a), each = b)

      sums <- rowsum(z, group = batch_id, reorder = FALSE)
      sums / b
    })
  )

  centred <- sweep(batch_means, 2, global_mean, FUN = "-")

  b * crossprod(centred) / (m * a - 1)
}
