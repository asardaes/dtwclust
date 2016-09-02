#' Keogh's DTW lower bound
#'
#' This function calculates a lower bound (LB) on the Dynamic Time Warp (DTW) distance between two time
#' series. It uses a Sakoe-Chiba constraint.
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
#' Keogh E and Ratanamahatana CA (2005). ``Exact indexing of dynamic time warping.'' \emph{Knowledge
#' and information systems}, \strong{7}(3), pp. 358-386.
#'
#' @examples
#'
#' # Sample data
#' data(uciCT)
#'
#' # Lower bound distance between two series
#' d.lbk <- lb_keogh(CharTraj[[1]], CharTraj[[2]], window.size = 20)$d
#'
#' # Corresponding true DTW distance
#' d.dtw <- dtw(CharTraj[[1]], CharTraj[[2]],
#'              window.type = "slantedband", window.size = 20)$distance
#'
#' d.lbk <= d.dtw
#'
#' # Calculating the LB between several time series using the 'proxy' package
#' # (notice how both argments must be lists)
#' D.lbk <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "LB_Keogh",
#'                      window.size = 20, norm = "L2")
#'
#' # Corresponding true DTW distance
#' # (see dtwclust-package description for an explanation of DTW2)
#' D.dtw <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "DTW2",
#'                      window.type = "slantedband", window.size = 20)
#'
#' D.lbk <= D.dtw
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
#' @return A list with: \itemize{
#'   \item \code{d}: The lower bound of the DTW distance.
#'   \item \code{upper.env}: The time series of \code{y}'s upper envelope.
#'   \item \code{lower.env}: The time series of \code{y}'s lower envelope.
#' }
#'
#' @export
#'

lb_keogh <- function(x, y, window.size = NULL, norm = "L1",
                     lower.env = NULL, upper.env = NULL, force.symmetry = FALSE) {

     norm <- match.arg(norm, c("L1", "L2"))

     consistency_check(x, "ts")

     if (is.null(lower.env) || is.null(upper.env)) {
          consistency_check(y, "ts")

          if (length(x) != length(y))
               stop("The series must have the same length")

          window.size <- consistency_check(window.size, "window")

          if (window.size > length(x))
               stop("The width of the window should not exceed the length of the series")
     }

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

     D <- rep(0, length(x))

     ind1 <- x > upper.env
     D[ind1] <- x[ind1] - upper.env[ind1]
     ind2 <- x < lower.env
     D[ind2] <- lower.env[ind2] - x[ind2]

     d <- switch(EXPR = norm,
                 L1 = sum(D),
                 L2 = sqrt(sum(D^2)))

     if (force.symmetry) {
          d2 <- lb_keogh(x = y, y = x, window.size = window.size, norm = norm)

          if (d2$d > d) {
               d <- d2$d
               lower.env = d2$lower.env
               upper.env = d2$upper.env
          }
     }

     ## Finish
     list(d = d,
          upper.env = upper.env,
          lower.env = lower.env)
}

# ========================================================================================================
# Loop without using native 'proxy' looping (to avoid multiple calculations of the envelope)
# ========================================================================================================

lb_keogh_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1",
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
          lower.env <- split_parallel(lower.env)
          upper.env <- split_parallel(upper.env)
     }

     if (pairwise) {
          D <- foreach(x = x, lower.env = lower.env, upper.env = upper.env,
                       .combine = c,
                       .multicombine = TRUE) %dopar% {
                            mapply(upper.env, lower.env, x,
                                   FUN = function(u, l, x) {
                                        lb_keogh(x,
                                                 norm = norm,
                                                 lower.env = l,
                                                 upper.env = u)$d
                                   })
                       }

          attr(D, "class") <- "pairdist"

     } else {
          D <- foreach(x = x,
                       .combine = rbind,
                       .multicombine = TRUE) %dopar% {
                            ret <- lapply(X = x, U = upper.env, L = lower.env,
                                          FUN = function(x, U, L) {
                                               ## This will return one row of the distance matrix
                                               D <- mapply(U, L,
                                                           MoreArgs = list(x = x),
                                                           FUN = function(u, l, x) {
                                                                lb_keogh(x,
                                                                         norm = norm,
                                                                         lower.env = l,
                                                                         upper.env = u)$d
                                                           })
                                               D
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

     attr(D, "method") <- "LB_Keogh"

     D
}
