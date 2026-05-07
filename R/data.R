#' Network tomography example: A6 around Leicester
#'
#' Data that was used in Hazelton et al. (2021) about dymanical lattice basis samplers
#'
#' @format ## `highway`
#' A list containing all that is needed to run the examples of the paper:
#' \describe{
#'   \item{A}{Configuration Matrix}
#'   \item{y}{Observation vector}
#'   \item{lambda}{Traffic rates for the different routes through A6}
#'   \item{LatticeBasis}{Lattice basis found using the column partition k^(-1)(1) = 1,2,3,4,5,6,7}
#'   \item{MarkovBasis}{A Markov basis that contains the lattice basis used for the extension samplers}
#' }
#' @source <doi: 10.1093/biomet/asaa083>
"highway"

#' 3-way contingency table: Book-crossing data
#'
#' Data that was used in Hazelton et al. (2021) and collected by Ziegler et al. (2005)
#'
#' @format ## `contingencyTable`
#' A list containing all that is needed to run the examples of the paper:
#' \describe{
#'   \item{A}{Configuration Matrix}
#'   \item{y}{Observation vector}
#'   \item{lambda}{Poisson rates for the cells of the table}
#'   \item{LatticeBasis}{Lattice basis found using the column partition as discussed in Extension Sampler chapter}
#'   \item{BasicMoves}{Collection of moves that forms a 1-Markov basis for this example}
#'   \item{MarkovBasis}{A Markov basis that contains the lattice basis and the collection of basic moves.}
#' }
#' @source <doi: 10.1093/biomet/asaa083>
"contingencyTable"
