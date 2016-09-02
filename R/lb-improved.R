#' Lemire's improved DTW lower bound
#'
#' This function calculates an improved lower bound (LB) on the Dynamic Time Warp (DTW) distance between two
#' time series. It uses a Sakoe-Chiba constraint.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10(2) + 1 = 21} observations falling within
#' the window.
#'
#' The reference time series should go in \code{x}, whereas the query time series should go in \code{y}.
#'
#' @note
#'
#' The lower bound is defined for time series of equal length only and is \strong{not} symmetric.
#'
#' If you wish to calculate the lower bound between several time series, it would be better to use the version
#' registered with the \code{proxy} package, since it includes some small optimizations. The convention
#' mentioned above for references and queries still holds. See the examples.
#'
#' @references
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .''
#' \emph{Pattern Recognition}, \strong{42}(9), pp. 2169 - 2180. ISSN 0031-3203,
#' \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030},
#' \url{http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
#'
#' @examples
#'
#' # Sample data
#' data(uciCT)
#'
#' # Lower bound distance between two series
#' d.lbi <- lb_improved(CharTraj[[1]], CharTraj[[2]], window.size = 20)
#'
#' # Corresponding true DTW distance
#' d.dtw <- dtw(CharTraj[[1]], CharTraj[[2]],
#'              window.type = "slantedband", window.size = 20)$distance
#'
#' d.lbi <= d.dtw
#'
#' # Calculating the LB between several time series using the 'proxy' package
#' # (notice how both argments must be lists)
#' D.lbi <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "LB_Improved",
#'                      window.size = 20, norm = "L2")
#'
#' # Corresponding true DTW distance
#' # (see dtwclust documentation for an explanation of DTW2)
#' D.dtw <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "DTW2",
#'                      window.type = "slantedband", window.size = 20)
#'
#' D.lbi <= D.dtw
#'
#' @param x A time series (reference).
#' @param y A time series with the same length as \code{x} (query).
#' @param window.size Window size for envelope calculation. See details.
#' @param norm Vector norm. Either \code{"L1"} for Manhattan distance or \code{"L2"} for Euclidean.
#' @param lower.env Optionally, a pre-computed lower envelope for \strong{\code{y}} can be provided
#' (non-proxy version only).
#' @param upper.env Optionally, a pre-computed upper envelope for \strong{\code{y}} can be provided
#' (non-proxy version only).
#' @param force.symmetry If \code{TRUE}, a second lower bound is calculated by swapping \code{x} and
#' \code{y}, and whichever result has a \emph{higher} distance value is returned. The proxy version
#' can only work if a square matrix is obtained, but use carefully.
#'
#' @return The improved lower bound for the DTW distance.
#'
#' @export

lb_improved <- function(x, y, window.size = NULL, norm = "L1",
                        lower.env = NULL, upper.env = NULL, force.symmetry = FALSE) {

     norm <- match.arg(norm, c("L1", "L2"))

     consistency_check(x, "ts")
     consistency_check(y, "ts")

     if (length(x) != length(y))
          stop("The series must have the same length")

     window.size <- consistency_check(window.size, "window")

     if (window.size > length(x))
          stop("The width of the window should not exceed the length of the series")

     ## LB Keogh first

     ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runminmax'
     if (is.null(lower.env) && is.null(upper.env)) {
          envelopes <- call_envelop(y, window.size*2L + 1L)
          lower.env <- envelopes$min
          upper.env <- envelopes$max

     } else if (is.null(lower.env)) {
          lower.env <- caTools::runmin(y, window.size*2L + 1L)

     } else if (is.null(upper.env)) {
          upper.env <- caTools::runmax(y, window.size*2L + 1L)
     }

     if (length(lower.env) != length(x))
          stop("Length mismatch between 'x' and the lower envelope")

     if (length(upper.env) != length(x))
          stop("Length mismatch between 'x' and the upper envelope")

     ind1 <- x > upper.env
     ind2 <- x < lower.env

     H <- x
     H[ind1] <- upper.env[ind1]
     H[ind2] <- lower.env[ind2]

     d1 <- abs(x - H)

     ## From here on is Lemire's improvement
     EH <- call_envelop(H, window.size*2+1)

     ind3 <- y > EH$max
     ind4 <- y < EH$min

     H2 <- y
     H2[ind3] <- EH$max[ind3]
     H2[ind4] <- EH$min[ind4]

     d2 <- abs(y - H2)

     ## LB_Improved is defined as root-p of the sum of LB_Keoghs^p
     d <- switch(EXPR = norm,
                 L1 = sum(d1) + sum(d2),
                 L2 = sqrt(sum(d1^2) + sum(d2^2))
     )

     if (force.symmetry) {
          d2 <- lb_improved(x = y, y = x, window.size = window.size, norm = norm)

          if (d2 > d) d <- d2
     }

     ## Finish
     d
}

# ========================================================================================================
# Loop without using native 'proxy' looping (to avoid multiple calculations of the envelope)
# ========================================================================================================

lb_improved_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1",
                              force.symmetry = FALSE, pairwise = FALSE, error.check = TRUE) {

     norm <- match.arg(norm, c("L1", "L2"))

     if (error.check)
          window.size <- consistency_check(window.size, "window")

     x <- consistency_check(x, "tsmat")

     if (error.check)
          consistency_check(x, "tslist")

     if (window.size > length(x[[1L]]))
          stop("Window size should not exceed length of the time series")

     if (is.null(y)) {
          y <- x

     } else {
          y <- consistency_check(y, "tsmat")

          if (error.check)
               consistency_check(y, "tslist")

          if (window.size > length(y[[1L]]))
               stop("Window size should not exceed length of the time series")
     }

     ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runminmax'
     envelops <- lapply(y, function(s) { call_envelop(s, window.size*2L + 1L) })

     lower.env <- lapply(envelops, "[[", "min")
     upper.env <- lapply(envelops, "[[", "max")

     check_parallel()

     x <- split_parallel(x)

     if (pairwise) {
          y <- split_parallel(y)
          lower.env <- split_parallel(lower.env)
          upper.env <- split_parallel(upper.env)
     }

     if (pairwise) {
          D <- foreach(x = x, y = y, lower.env = lower.env, upper.env = upper.env,
                       .combine = c,
                       .multicombine = TRUE,
                       .packages = "dtwclust") %dopar% {
                            mapply(upper.env, lower.env, y, x,
                                   FUN = function(u, l, y, x) {
                                        lb_improved(x, y,
                                                    window.size = window.size,
                                                    norm = norm,
                                                    lower.env = l,
                                                    upper.env = u)
                                   })
                       }

          attr(D, "class") <- "pairdist"

     } else {
          D <- foreach(x = x,
                       .combine = rbind,
                       .multicombine = TRUE,
                       .packages = "dtwclust") %dopar% {
                            ret <- lapply(X = x, U = upper.env, L = lower.env, Y = y,
                                          FUN = function(x, U, L, Y) {
                                               ## This will return one row of the distance matrix
                                               mapply(U, L, Y,
                                                      MoreArgs = list(x = x),
                                                      FUN = function(u, l, y, x) {
                                                           lb_improved(x, y,
                                                                       window.size = window.size,
                                                                       norm = norm,
                                                                       lower.env = l,
                                                                       upper.env = u)
                                                      })
                                          })

                            do.call(rbind, ret)
                       }

          attr(D, "class") <- "crossdist"
     }

     if (force.symmetry && !pairwise) {
          if (nrow(D) != ncol(D)) {
               warning("Unable to force symmetry. Resulting distance matrix is not square.")

          } else {
               ind.tri <- lower.tri(D)

               new.low.tri.vals <- t(D)[ind.tri]
               indCorrect <- D[ind.tri] > new.low.tri.vals
               new.low.tri.vals[indCorrect] <- D[ind.tri][indCorrect]

               D[ind.tri] <- new.low.tri.vals
               D <- t(D)
               D[ind.tri] <- new.low.tri.vals
          }
     }

     attr(D, "method") <- "LB_Improved"

     D
}
