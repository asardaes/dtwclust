#' Shape-based distance
#'
#' Distance based on coefficient-normalized cross-correlation as proposed by Papparizos and Gravano, 2015,
#' for the k-Shape clustering algorithm.
#'
#' This function works best if the series are \emph{z-normalized}. If not, at least they should have
#' corresponding amplitudes, since the values of the signal \strong{do} affect the outcome.
#'
#' If \code{x} and \code{y} do \strong{not} have the same length, it would be best if the longer sequence is
#' provided in \code{y}, because it will be shifted to match \code{x}. Anything before the matching point is
#' discarded and the series is padded with trailing zeros as needed.
#'
#' The output values lie between 0 and 2, with 0 indicating perfect similarity.
#'
#' @examples
#'
#' # load data
#' data(uciCT)
#'
#' # distance between series of different lengths
#' sbd <- SBD(CharTraj[[1]], CharTraj[[100]], znorm = TRUE)$dist
#'
#' # cross-distance matrix for series subset (notice the two-list input)
#' sbD <- proxy::dist(CharTraj[1:10], CharTraj[1:10], method = "SBD", znorm = TRUE)
#'
#' @seealso
#'
#' \code{\link{NCCc}}, \code{\link{shape_extraction}}
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.'' In \emph{Proceedings of the 2015
#' ACM SIGMOD International Conference on Management of Data}, series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{
#' http://doi.org/10.1145/2723372.2737793}.
#'
#' @param x A time series.
#' @param y Another time series.
#' @param znorm Should each series be z-normalized before calculating the distance?
#'
#' @return A list with: \itemize{
#'   \item \code{dist}: The distance between \code{x} and \code{y}.
#'   \item \code{yshift}: A shifted version of \code{y} so that it optimally mathces \code{x}.
#' }
#'
#' @export

SBD <- function(x, y, znorm = FALSE) {

     x <- as.numeric(x)
     y <- as.numeric(y)

     consistency_check(x, "ts")
     consistency_check(y, "ts")

     if (znorm)
          CCseq <- NCCc(zscore(x), zscore(y))
     else
          CCseq <- NCCc(x,y)

     m <- max(CCseq)
     d <- which.max(CCseq)
     n <- length(y)

     shift <- d - max(length(x), length(y))

     if (shift < 0)
          yshift <- c( y[(-shift+1):n], rep(0, -shift) )
     else
          yshift <- c( rep(0, shift), y[1:(n-shift)] )

     dist <- 1 - m

     list(dist = dist,
          yshift = yshift)
}

# ========================================================================================================
# Wrapper for proxy::dist
# ========================================================================================================

SBD.proxy <- function(x, y, znorm = FALSE) {

     x <- as.numeric(x)
     y <- as.numeric(y)

     consistency_check(x, "ts")
     consistency_check(y, "ts")

     if (znorm)
          CCseq <- NCCc(zscore(x), zscore(y))
     else
          CCseq <- NCCc(x,y)

     dist <- 1 - max(CCseq)

     dist
}
