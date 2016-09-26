#' Shape-based distance
#'
#' Distance based on coefficient-normalized cross-correlation as proposed by Paparrizos and Gravano, 2015,
#' for the k-Shape clustering algorithm.
#'
#' This distance works best if the series are \emph{z-normalized}. If not, at least they should have
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
#' This distance is calculated with help of the Fast Fourier Transform, so it can be sensitive to numerical
#' precision. Results could vary slightly between 32 and 64 bit architectures.
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
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In \emph{Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data},
#' series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' @param x,y A time series.
#' @param znorm Logical. Should each series be z-normalized before calculating the distance?
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
     shift <- which.max(CCseq) - max(nx, ny)

     if (is.null(flip)) {
          if (shift < 0L)
               yshift <- y[(-shift + 1L):ny]
          else
               yshift <- c( rep(0, shift), y )

     } else {
          ## Remember, if I flipped them, then I have to shift what is now saved in 'x'
          if (shift < 0L)
               yshift <- c( rep(0, -shift), x )
          else
               yshift <- x[(shift + 1L):ny]
     }

     nys <- length(yshift)

     if (nys < nx)
          yshift <- c( yshift, rep(0, nx-nys) )
     else
          yshift <- yshift[1L:nx]

     list(dist = 1 - m, yshift = yshift)
}

# ========================================================================================================
# Wrapper for proxy::dist
# ========================================================================================================

SBD.proxy <- function(x, y = NULL, znorm = FALSE, error.check = TRUE, pairwise = FALSE, ...) {
     x <- consistency_check(x, "tsmat")

     if (error.check)
          consistency_check(x, "vltslist")

     if(znorm) x <- zscore(x)

     if (is.null(y)) {
          y <- x

     } else {
          y <- consistency_check(y, "tsmat")

          if (error.check)
               consistency_check(y, "vltslist")

          if(znorm) y <- zscore(y)
     }

     retclass <- "crossdist"

     ## Precompute FFTs, padding with zeros as necessary, which will be compensated later
     L <- max(lengths(x)) + max(lengths(y)) - 1L
     fftlen <- stats::nextn(L, 2L)

     fftx <- lapply(x, function(u) { stats::fft(c(u, rep(0, fftlen - length(u)))) })

     ffty <- lapply(y, function(v) { stats::fft(c(v, rep(0, fftlen - length(v)))) })

     ## Register doSEQ if necessary
     check_parallel()

     x <- split_parallel(x)
     fftx <- split_parallel(fftx)

     ## Calculate distance matrix
     if (pairwise) {
          y <- split_parallel(y)
          ffty <- split_parallel(ffty)

          if (length(lengths(x)) != length(lengths(y)) || lengths(x) != lengths(y))
               stop("Pairwise distances require the same amount of series in 'x' and 'y'")

          D <- foreach(x = x, fftx = fftx, y = y, ffty = ffty,
                       .combine = c,
                       .multicombine = TRUE,
                       .export = "lnorm",
                       .packages = "stats") %dopar% {
                            mapply(y, ffty, x, fftx,
                                   FUN = function(y, ffty, x, fftx) {

                                        ## Manually normalize by length
                                        CCseq <- Re(stats::fft(fftx * Conj(ffty), inverse = TRUE)) / length(fftx)

                                        ## Truncate to correct length
                                        CCseq <- c(CCseq[(length(ffty) - length(y) + 2L):length(CCseq)],
                                                   CCseq[1L:length(x)])

                                        CCseq <- CCseq / (lnorm(x, 2) * lnorm(y, 2))

                                        dd <- 1 - max(CCseq)

                                        dd
                                   })
                       }

          retclass <- "pairdist"

     } else {
          D <- foreach(x = x, fftx = fftx,
                       .combine = rbind,
                       .multicombine = TRUE,
                       .export = "lnorm",
                       .packages = "stats") %dopar% {
                            ret <- mapply(x, fftx,
                                          MoreArgs = list(Y = y, FFTY = ffty),
                                          SIMPLIFY = FALSE,
                                          FUN = function(x, fftx, Y, FFTY) {
                                               mapply(Y, FFTY,
                                                      MoreArgs = list(x = x, fftx = fftx),
                                                      FUN = function(y, ffty, x, fftx) {
                                                           ## Manually normalize by length
                                                           CCseq <- Re(stats::fft(fftx * Conj(ffty),
                                                                                  inverse = TRUE)) / length(fftx)

                                                           ## Truncate to correct length
                                                           CCseq <- c(CCseq[(length(ffty) - length(y) + 2L):length(CCseq)],
                                                                      CCseq[1L:length(x)])

                                                           CCseq <- CCseq / (lnorm(x, 2) * lnorm(y, 2))

                                                           dd <- 1 - max(CCseq)

                                                           dd
                                                      })
                                          })

                            do.call(rbind, ret)
                       }
     }

     class(D) <- retclass
     attr(D, "method") <- "SBD"

     D
}
