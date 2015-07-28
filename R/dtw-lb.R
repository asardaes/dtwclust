#' DTW calculation guided by Lemire's lower bound
#'
#' Calculation of a distance matrix with the Dynamic Time Warp (DTW) distance guided by Lemire's lower bound
#' (LB).
#'
#' This function first calculates an initial estimate of a distance matrix between two sets of time series
#' using Lemire's improved lower bound. Afterwards, it uses the estimate to calculate the true DTW distances
#' of \emph{only} the nearest neighbors for each series in \code{x}.
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
#' The function takes care of this discrepancy if needed, but you should be careful if you are
#' testing things manually.
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
     if (class(x) == "matrix")
          X <- lapply(seq_len(nrow(x)), function(i) x[i,])
     else
          X <- x

     if (!is.null(y) && class(y) == "matrix")
          Y <- lapply(seq_len(nrow(y)), function(i) y[i,])
     else
          Y <- X

     ## Initial estimate
     d <- lb_improved_loop(X, Y, window.size = window.size*2+1, norm = norm, error.check = error.check)

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
                              window.type = "sakoechiba", window.size = window.size),

                         L2 = proxy::dist(X[indNew], Y[indNN[indNew]], pairwise = TRUE,
                                          method = "DTW2",
                                          window.type = "sakoechiba", window.size = window.size))

          indD <- indNN + singleIndexing
          d[indD[indNew]] <- dSub

          new.indNN <- apply(d, 2, which.min)
     }

     ## Transpose again for final result
     attr(d, "method") <- "DTW_LB"
     return(t(d))
}
