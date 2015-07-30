#' DTW calculation guided by Lemire's lower bound (LB_Improved)
#'
#' Calculation of a distance matrix with the Dynamic Time Warping (DTW) distance guided by Lemire's lower bound
#' (LB).
#'
#' This function first calculates an initial estimate of a distance matrix between two sets of time series
#' using Lemire's improved lower bound. Afterwards, it uses the estimate to calculate the true DTW distances
#' between \emph{only} the nearest neighbors of each series in \code{x} found in \code{y}. If only \code{x}
#' is provided, the distance matrix is calculated between all its time series. This could be useful in case
#' one is interested in only the nearest neighbor of one or more series among a dataset.
#'
#' Because of the way the different functions being used here are implemented, there is a subtle but critical
#' mismatch in the way the window size is defined for DTW and the LB. The DTW calculation with
#' \code{\link[dtw]{dtw}} expects an \emph{even} \code{window.size} that represents the distance between the
#' diagonal and one of the edges of the window. The LB calculation expects an \emph{odd} window.size that
#' represents the whole window width to be used in the running max and min. The outcome of said running
#' functions are centered with respect to the window width.
#'
#' Therefore, if, for example, the DTW is calculated with a window of 10, the corresponding LB should
#' be calculated with \code{2*10 + 1 = 21}.
#'
#' This function expects the \code{window.size} for \code{DTW} and takes care of this discrepancy automatically.
#'
#' @note
#'
#' This function uses a lower bound that is only defined for time series of equal lengths.
#'
#' @seealso
#'
#' \code{\link{lb_improved}}
#'
#' @references
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .'' \emph{Pattern Recognition}, \strong{42}(9), pp.
#' 2169 - 2180. ISSN 0031-3203, \url{http://doi.org/http://dx.doi.org/10.1016/j.patcog.2008.11.030}, \url{
#' http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
#'
#' @examples
#'
#' # Load data
#' data(uciCT)
#'
#' # Reinterpolate to same length
#' data <- lapply(CharTraj, reinterpolate, newLength = 205)
#'
#' # Calculate the DTW distance between a certain subset aided with the lower bound
#' system.time(d <- dtw_lb(data[1:5], data[6:100], window.size = 20))
#'
#' # Nearest neighbors
#' NN1 <- apply(d, 1, which.min)
#'
#' # Calculate the DTW distances between all elements (about seven times slower)
#' system.time(d2 <- proxy::dist(data[1:5], data[6:100], method = "DTW",
#'                               window.type = "slantedband", window.size = 20))
#'
#' # Nearest neighbors
#' NN2 <- apply(d2, 1, which.min)
#'
#' all(NN1 == NN2)
#'
#' @param x A matrix where rows are time series, or a list of time series.
#' @param y An object similar to \code{x}.
#' @param window.size The window size to use with the DTW calculation. \strong{See details}.
#' @param norm Pointwise distance. Either \code{L1} for Manhattan distance or \code{L2} for Euclidean.
#' @param error.check Should inconsistencies in the data be checked?
#'
#' @return The distance matrix with class \code{crossdist}.
#'
#' @export
#' @importFrom dtw dtw

dtw_lb <- function(x, y = NULL, window.size = NULL, norm = "L1", error.check = TRUE) {

     if (is.null(window.size)) {
          stop("Please provide the 'window.size' parameter")
     }
     if (window.size%%2 != 0) {
          stop("For the Sakoe-Chiba band, the window must be symmetric and window.size must be even")
     }
     if (window.size <= 1) {
          stop("Window width must be larger than 1")
     }

     norm <- match.arg(norm, c("L1", "L2"))

     ## For looping convenience
     if (is.matrix(x))
          X <- lapply(seq_len(nrow(x)), function(i) x[i,])
     else if (is.numeric(x))
          X <- list(x)
     else if (is.list(x))
          X <- x
     else
          stop("Unsupported type for x")

     if (!is.null(y)) {
          if (is.matrix(y))
               Y <- lapply(seq_len(nrow(y)), function(i) y[i,])
          else if (is.numeric(y))
               Y <- list(y)
          else if (is.list(y))
               Y <- y
          else
               stop("Unsupported type for y")

     } else {
          Y <- X
     }

     ## Initial estimate
     d <- proxy::dist(X, Y, method = "LBI",
                      window.size = window.size*2+1,
                      norm = norm,
                      error.check = error.check)

     ## For indexing convenience
     d <- t(d)
     singleIndexing <- seq(from=0, by=nrow(d), length.out=ncol(d))

     ## Update with DTW
     ## NOTE: the 'window.size' definition varies with respect to the LB, see above

     new.indNN <- apply(d, 2, which.min) # index of nearest neighbors
     indNN <- new.indNN + 1

     while (any(new.indNN != indNN)) {
          indNew <- new.indNN != indNN
          indNN <- new.indNN

          dSub <- switch(EXPR = norm,

                         L1 = proxy::dist(X[indNew], Y[indNN[indNew]], pairwise = TRUE,
                                          method = "DTW", dist.method = "L1",
                                          window.type = "slantedband", window.size = window.size),

                         L2 = proxy::dist(X[indNew], Y[indNN[indNew]], pairwise = TRUE,
                                          method = "DTW2",
                                          window.type = "slantedband", window.size = window.size))

          indD <- indNN + singleIndexing
          d[indD[indNew]] <- dSub

          new.indNN <- apply(d, 2, which.min)
     }

     ## Transpose again for final result
     attr(d, "method") <- "DTW_LB"
     t(d)
}
