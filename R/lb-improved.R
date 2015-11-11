#' Lemire's improved DTW lower bound
#'
#' This function calculates a lower bound (LB) on the Dynamic Time Warp (DTW) distance between two time
#' series. It uses a Sakoe-Chiba constraint.
#'
#' The lower bound is defined for time series of equal length only.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10*2 + 1 = 21} observations falling within
#' the window.
#'
#' @note
#'
#' If you wish to calculate the lower bound between several time series, it would be better to use the version
#' registered with the 'proxy' package, since it includes some small optimizations. See the examples.
#'
#' @references
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .'' \emph{Pattern Recognition}, \strong{42}(9), pp.
#' 2169 - 2180. ISSN 0031-3203, \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030}, \url{
#' http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
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
#' # (see dtwclust-package description for an explanation of DTW2)
#' D.dtw <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "DTW2",
#'                      window.type = "slantedband", window.size = 20)
#'
#' D.lbi <= D.dtw
#'
#' @param x A time series.
#' @param y A time series with the same length as \code{x}.
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

     if (length(x) != length(y)) {
          stop("The series must have the same length")
     }

     window.size <- consistency_check(window.size, "window")

     if (window.size > length(x))
          stop("The width of the window should not exceed the length of the series")

     ## LB Keogh first

     ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runmax/min'
     if (is.null(lower.env)) {
          lower.env <- caTools::runmin(y, window.size*2+1, endrule="constant")
     } else {
          if (length(lower.env) != length(y))
               stop("Length mismatch between 'y' and its lower envelope")
     }

     if (is.null(upper.env)) {
          upper.env <- caTools::runmax(y, window.size*2+1, endrule="constant")
     } else {
          if (length(upper.env) != length(y))
               stop("Length mismatch between 'y' and its upper envelope")
     }

     ind1 <- which(x > upper.env)
     ind2 <- which(x < lower.env)

     H <- x
     H[ind1] <- upper.env[ind1]
     H[ind2] <- lower.env[ind2]

     d1 <- abs(x-H)

     ## From here on is Lemire's improvement
     UH <- caTools::runmax(H, window.size*2+1, endrule="constant")
     LH <- caTools::runmin(H, window.size*2+1, endrule="constant")

     ind3 <- which(y > UH)
     ind4 <- which(y < LH)

     H2 <- y
     H2[ind3] <- UH[ind3]
     H2[ind4] <- LH[ind4]

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

lb_improved_loop <- function(x, y = NULL, ...) {

     ARGS <- list(...)
     window.size <- ARGS$window.size
     error.check <- ARGS$error.check
     force.symmetry <- ARGS$force.symmetry
     norm <- ARGS$norm

     if (is.null(error.check))
          error.check <- TRUE
     if (is.null(force.symmetry))
          force.symmetry <- FALSE

     if (is.null(norm))
          norm <- "L1"
     else
          norm <- match.arg(norm, c("L1", "L2"))

     if (error.check)
          window.size <- consistency_check(window.size, "window")

     x <- consistency_check(x, "tsmat")

     if (error.check)
          consistency_check(x, "tslist")

     if (window.size > length(x[[1]]))
          stop("Window size should not exceed length of the time series")

     if (is.null(y)) {
          y <- x

     } else {

          y <- consistency_check(y, "tsmat")

          if (error.check)
               consistency_check(y, "tslist")

          if (window.size > length(y[[1]]))
               stop("Window size should not exceed length of the time series")
     }

     ## from 'caTools' package
     ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runmax/min'
     upper.env <- lapply(y, runmax, k=window.size*2+1, endrule="constant")
     lower.env <- lapply(y, runmin, k=window.size*2+1, endrule="constant")

     DD <- sapply(X=x, U=upper.env, L=lower.env, Y=y,
                  FUN = function(x, ...) {
                       U <- list(...)$U
                       L <- list(...)$L
                       Y <- list(...)$Y

                       ## This will return one column of the distance matrix
                       D <- mapply(U, L, Y, MoreArgs=list(x=x),
                                   FUN = function(u, l, y, x) {

                                        ## LB Keogh
                                        ind1 <- which(x > u)
                                        ind2 <- which(x < l)

                                        H <- x
                                        H[ind1] <- u[ind1]
                                        H[ind2] <- l[ind2]

                                        d1 <- abs(x-H)

                                        ## Lemire's improvement
                                        UH <- caTools::runmax(H, window.size*2+1, endrule="constant")
                                        LH <- caTools::runmin(H, window.size*2+1, endrule="constant")

                                        ind3 <- which(y > UH)
                                        ind4 <- which(y < LH)

                                        H2 <- y
                                        H2[ind3] <- UH[ind3]
                                        H2[ind4] <- LH[ind4]

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

     if (force.symmetry) {
          if (nrow(DD) != ncol(DD)) {
               warning("Unable to force symmetry. Resulting distance matrix is not square.")
          } else {
               ind.tri <- lower.tri(DD)

               new.low.tri.vals <- t(DD)[ind.tri]
               indCorrect <- DD[ind.tri] > new.low.tri.vals
               new.low.tri.vals[indCorrect] <- DD[ind.tri][indCorrect]

               DD[ind.tri] <- new.low.tri.vals
               DD <- t(DD)
               DD[ind.tri] <- new.low.tri.vals
          }
     }

     attr(DD, "class") <- "crossdist"
     attr(DD, "method") <- "LB_Improved"

     t(DD)
}
