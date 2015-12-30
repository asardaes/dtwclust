#' Shape-based distance
#'
#' Distance based on coefficient-normalized cross-correlation as proposed by Papparizos and Gravano, 2015,
#' for the k-Shape clustering algorithm.
#'
#' This function works best if the series are \emph{z-normalized}. If not, at least they should have
#' corresponding amplitudes, since the values of the signals \strong{do} affect the outcome.
#'
#' If \code{x} and \code{y} do \strong{not} have the same length, it would be best if the longer sequence is
#' provided in \code{y}, because it will be shifted to match \code{x}. Anything before the matching point is
#' discarded and the series is padded with trailing zeros as needed.
#'
#' The output values lie between 0 and 2, with 0 indicating perfect similarity.
#'
#' @note
#'
#' If you wish to calculate the distance between several time series, it would be better to use the version
#' registered with the \code{proxy} package, since it includes some small optimizations. See the examples.
#'
#' However, because of said optimizations and the way \code{proxy}'s \code{\link[proxy]{dist}} works, the
#' latter's \code{pairwise} argument will not work with this distance. You can use the custom argument
#' \code{force.pairwise} to get the correct result.
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
#' @param x,y A time series.
#' @param znorm Should each series be z-normalized before calculating the distance?
#'
#' @return A list with: \itemize{
#'   \item \code{dist}: The shape-based distance between \code{x} and \code{y}.
#'   \item \code{yshift}: A shifted version of \code{y} so that it optimally matches \code{x}.
#' }
#'
#' @export
#'

SBD <- function(x, y, znorm = FALSE) {

     x <- as.numeric(x)
     y <- as.numeric(y)

     consistency_check(x, "ts")
     consistency_check(y, "ts")

     nx <- length(x)
     ny <- length(y)

     if (nx > ny) {
          ## The order in which I provide the arguments to NCCc affects 'shift'
          flip <- x
          x <- y
          y <- flip

     } else {
          flip <- NULL
     }

     if (znorm)
          CCseq <- NCCc(zscore(x), zscore(y))
     else
          CCseq <- NCCc(x,y)

     m <- max(CCseq)
     d <- which.max(CCseq)

     shift <- d - max(nx, ny)

     if (is.null(flip)) {
          if (shift < 0)
               yshift <- y[(-shift+1):ny]
          else
               yshift <- c( rep(0, shift), y )

     } else {
          ## Remember, if I flipped them, then I have to shift what is now saved in 'x'
          if (shift < 0)
               yshift <- c( rep(0, -shift), x )
          else
               yshift <- x[(shift+1):ny]
     }

     nys <- length(yshift)

     if (nys < nx)
          yshift <- c( yshift, rep(0, nx-nys) )
     else
          yshift <- yshift[1:nx]

     dist <- 1 - m

     list(dist = dist,
          yshift = yshift)
}

# ========================================================================================================
# Wrapper for proxy::dist
# ========================================================================================================

SBD.proxy <- function(x, y = NULL, znorm = FALSE, error.check = TRUE, force.pairwise = FALSE) {

     x <- consistency_check(x, "tsmat")

     if (error.check)
          consistency_check(x, "vltslist")

     if (is.null(y)) {
          if(znorm) {
               x <- lapply(x, function(xx) as.numeric(zscore(xx)))
          }

          y <- x

     } else {
          y <- consistency_check(y, "tsmat")

          if (error.check)
               consistency_check(y, "vltslist")

          if(znorm) {
               x <- lapply(x, function(xx) as.numeric(zscore(xx)))
               y <- lapply(y, function(yy) as.numeric(zscore(yy)))
          }
     }

     ## Precompute FFTs, padding with zeros as necessary, which will be compensated later
     L <- max(sapply(x, length)) + max(sapply(y, length)) - 1
     fftlen <- stats::nextn(L, 2)

     fftx <- lapply(x, function(u) {
          stats::fft(c(u, rep(0L, fftlen-length(u))))
     })

     ffty <- lapply(y, function(v) {
          stats::fft(c(v, rep(0L, fftlen-length(v))))
     })

     check_parallel()

     x <- split_parallel(x)
     fftx <- split_parallel(fftx)

     if (force.pairwise) {
          y <- split_parallel(y)
          ffty <- split_parallel(ffty)
     }

     ## Calculate distance matrix
     if (force.pairwise) {
          D <- foreach(x = x, fftx = fftx, y = y, ffty = ffty,
                       .combine = c,
                       .multicombine = TRUE,
                       .packages = "stats") %dopar% {
                            mapply(y, ffty, x, fftx,
                                   FUN = function(y, ffty, x, fftx) {

                                        CCseq <- Re(stats::fft(fftx * Conj(ffty), inverse = TRUE)) / length(fftx)

                                        ## Truncate to correct length
                                        CCseq <- c(CCseq[(length(ffty)-length(y)+2):length(CCseq)],
                                                   CCseq[1:length(x)])
                                        CCseq <- CCseq / (sqrt(crossprod(x)) * sqrt(crossprod(y)))

                                        dd <- 1 - max(CCseq)

                                        dd
                                   })
                       }

          attr(D, "class") <- "pairdist"

     } else {
          D <- foreach(x = x, fftx = fftx,
                       .combine = cbind,
                       .multicombine = TRUE,
                       .packages = "stats") %dopar% {
                            mapply(x, fftx, MoreArgs = list(Y = y, FFTY = ffty),
                                   FUN = function(x, fftx, Y, FFTY) {

                                        d <- mapply(Y, FFTY, MoreArgs = list(x = x, fftx = fftx),
                                                    FUN = function(y, ffty, x, fftx) {

                                                         CCseq <- Re(stats::fft(fftx * Conj(ffty), inverse = TRUE)) / length(fftx)

                                                         ## Truncate to correct length
                                                         CCseq <- c(CCseq[(length(ffty)-length(y)+2):length(CCseq)],
                                                                    CCseq[1:length(x)])
                                                         CCseq <- CCseq / (sqrt(crossprod(x)) * sqrt(crossprod(y)))

                                                         dd <- 1 - max(CCseq)

                                                         dd
                                                    })

                                        d
                                   })
                       }

          attr(D, "class") <- "crossdist"
          D <- t(D)
     }

     attr(D, "method") <- "SBD"

     D
}
