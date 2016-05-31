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
#' registered with the \code{proxy} package, since it includes some small optimizations. See the examples.
#'
#' However, because of said optimizations and the way \code{proxy}'s \code{\link[proxy]{dist}} works, the
#' latter's \code{pairwise} argument will not work with this distance. You can use the custom argument
#' \code{force.pairwise} to get the correct result.
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
#' @param norm Pointwise distance. Either \code{L1} for Manhattan distance or \code{L2} for Euclidean.
#' @param lower.env Optionally, a pre-computed lower envelope for \strong{\code{y}} can be provided.
#' @param upper.env Optionally, a pre-computed upper envelope for \strong{\code{y}} can be provided.
#'
#' @return The improved lower bound for the DTW distance.
#'
#' @export

lb_improved <- function(x, y, window.size = NULL, norm = "L1", lower.env = NULL, upper.env = NULL) {

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

     if (length(lower.env) != length(y))
          stop("Length mismatch between 'y' and its lower envelope")

     if (length(upper.env) != length(y))
          stop("Length mismatch between 'y' and its upper envelope")

     ind1 <- x > upper.env
     ind2 <- x < lower.env

     H <- x
     H[ind1] <- upper.env[ind1]
     H[ind2] <- lower.env[ind2]

     d1 <- abs(x-H)

     ## From here on is Lemire's improvement
     EH <- call_envelop(H, window.size*2+1)

     ind3 <- y > EH$max
     ind4 <- y < EH$min

     H2 <- y
     H2[ind3] <- EH$max[ind3]
     H2[ind4] <- EH$min[ind4]

     d2 <- abs(y-H2)

     ## LB_Improved is defined as root-p of the sum of LB_Keoghs^p
     d <- switch(EXPR = norm,
                 L1 = sum(d1) + sum(d2),
                 L2 = sqrt(sum(d1^2) + sum(d2^2))
     )

     ## Finish
     d
}

# ========================================================================================================
# Loop without using native 'proxy' looping (to avoid multiple calculations of the envelope)
# ========================================================================================================

lb_improved_loop <- function(x, y = NULL, window.size = NULL, error.check = TRUE,
                             force.symmetry = FALSE, norm = "L1", force.pairwise = FALSE) {

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

     if (force.pairwise) {
          y <- split_parallel(y)
          lower.env <- split_parallel(lower.env)
          upper.env <- split_parallel(upper.env)
     }

     if (force.pairwise) {
          D <- foreach(x = x, y = y, lower.env = lower.env, upper.env = upper.env,
                       .combine = c,
                       .multicombine = TRUE,
                       .export = "call_envelop",
                       .packages = "dtwclust") %dopar% {
                            mapply(upper.env, lower.env, y, x,
                                   FUN = function(u, l, y, x) {

                                        ## LB Keogh
                                        ind1 <- x > u
                                        ind2 <- x < l

                                        H <- x
                                        H[ind1] <- u[ind1]
                                        H[ind2] <- l[ind2]

                                        d1 <- abs(x-H)

                                        ## Lemire's improvement
                                        EH <- call_envelop(H, window.size*2L + 1L)

                                        ind3 <- y > EH$max
                                        ind4 <- y < EH$min

                                        H2 <- y
                                        H2[ind3] <- EH$max[ind3]
                                        H2[ind4] <- EH$min[ind4]

                                        d2 <- abs(y-H2)

                                        ## calculate LB_Improved

                                        d <- switch(EXPR = norm,
                                                    L1 = sum(d1) + sum(d2),
                                                    L2 = sqrt(sum(d1^2) + sum(d2^2))
                                        )

                                        d
                                   })
                       }

          attr(D, "class") <- "pairdist"

     } else {
          D <- foreach(x = x,
                       .combine = cbind,
                       .multicombine = TRUE,
                       .export = "call_envelop",
                       .packages = "dtwclust") %dopar% {
                            sapply(X=x, U=upper.env, L=lower.env, Y=y,
                                   FUN = function(x, U, L, Y) {

                                        ## This will return one column of the distance matrix
                                        D <- mapply(U, L, Y, MoreArgs=list(x=x),
                                                    FUN = function(u, l, y, x) {

                                                         ## LB Keogh
                                                         ind1 <- x > u
                                                         ind2 <- x < l

                                                         H <- x
                                                         H[ind1] <- u[ind1]
                                                         H[ind2] <- l[ind2]

                                                         d1 <- abs(x-H)

                                                         ## Lemire's improvement
                                                         EH <- call_envelop(H, window.size*2L + 1L)

                                                         ind3 <- y > EH$max
                                                         ind4 <- y < EH$min

                                                         H2 <- y
                                                         H2[ind3] <- EH$max[ind3]
                                                         H2[ind4] <- EH$min[ind4]

                                                         d2 <- abs(y-H2)

                                                         ## calculate LB_Improved

                                                         d <- switch(EXPR = norm,
                                                                     L1 = sum(d1) + sum(d2),
                                                                     L2 = sqrt(sum(d1^2) + sum(d2^2))
                                                         )

                                                         d
                                                    })
                                        D
                                   })
                       }

          attr(D, "class") <- "crossdist"
          D <- t(D)
     }

     if (force.symmetry && !force.pairwise) {
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
