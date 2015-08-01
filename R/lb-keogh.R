#' Keogh's DTW lower bound
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
#' @references
#'
#' Keogh E and Ratanamahatana CA (2005). ``Exact indexing of dynamic time warping.'' \emph{Knowledge and information systems}, \strong{7}(3),
#' pp. 358-386.
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
#' d.dtw <- dtw(CharTraj[[1]], CharTraj[[2]], window.type = "slantedband", window.size = 20)$distance
#'
#' d.lbk <= d.dtw
#'
#' @param x A time series.
#' @param y A time series with the same length as \code{x}.
#' @param window.size Window size for envelope calculation. See details.
#' @param norm Pointwise distance. Either \code{L1} for Manhattan distance or \code{L2} for Euclidean.
#'
#' @return A list with: \itemize{
#'   \item \code{d}: The lower bound of the DTW distance.
#'   \item \code{upper.env}: The time series of the upper envelope.
#'   \item \code{lower.env}: The time series of the lower envelope.
#' }
#'
#' @export
#' @importFrom caTools runmax
#' @importFrom caTools runmin

lb_keogh <- function(x, y, window.size = NULL, norm = "L1") {

     norm <- match.arg(norm, c("L1", "L2"))

     consistency_check(x, "ts")
     consistency_check(y, "ts")

     if (length(x) != length(y)) {
          stop("The series must have the same length")
     }

     window.size <- consistency_check(window.size, "window")

     if (window.size > length(x))
          stop("The width of the window should not exceed the length of the series")

     # from 'caTools' package
     ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runmax/min'
     upper.env <- runmax(y, window.size*2+1, endrule="constant")
     lower.env <- runmin(y, window.size*2+1, endrule="constant")

     D <- rep(0, length(x))

     ind1 <- which(x > upper.env)
     D[ind1] <- x[ind1] - upper.env[ind1]
     ind2 <- which(x < lower.env)
     D[ind2] <- lower.env[ind2] - x[ind2]

     d <- switch(EXPR = norm,
                 L1 = sum(D),
                 L2 = sqrt(sum(D^2)))

     list(d = d,
          upper.env = upper.env,
          lower.env = lower.env)
}

# ========================================================================================================
# Loop without using native 'proxy' (to avoid multiple calculations of the envelope)
# - About 10 times faster
# ========================================================================================================

lb_keogh_loop <- function(x, y=NULL, ...) {

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

     ## For looping convenience
     if (is.matrix(x))
          x <- lapply(seq_len(nrow(x)), function(i) x[i,])
     else if (is.numeric(x))
          x <- list(x)
     else if (is.list(x))
          x <- x
     else
          stop("Unsupported type for x")

     if (error.check)
          consistency_check(x, "tslist")

     if (window.size > length(x[[1]]))
          stop("Window size should not exceed length of the time series")

     if (is.null(y)) {
          # from 'caTools' package
          ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runmax/min'
          upper.env <- lapply(x, runmax, k=window.size*2+1, endrule="constant")
          lower.env <- lapply(x, runmin, k=window.size*2+1, endrule="constant")

          DD <- sapply(X=x, U=upper.env, L=lower.env,
                       FUN = function(x, ...) {
                            U <- list(...)$U
                            L <- list(...)$L

                            ## This will return one column of the distance matrix
                            D <- mapply(U, L, MoreArgs=list(x=x),
                                        FUN = function(u, l, x) {

                                             D <- rep(0, length(x))

                                             ind1 <- which(x > u)
                                             D[ind1] <- x[ind1] - u[ind1]
                                             ind2 <- which(x < l)
                                             D[ind2] <- l[ind2] - x[ind2]

                                             d <- switch(EXPR = norm,
                                                         L1 = sum(D),
                                                         L2 = sqrt(sum(D^2)))

                                             d

                                        })

                            D
                       })

          if (force.symmetry) {
               ind.tri <- lower.tri(DD)

               new.low.tri.vals <- t(DD)[ind.tri]
               indCorrect <- DD[ind.tri] > new.low.tri.vals
               new.low.tri.vals[indCorrect] <- DD[ind.tri][indCorrect]

               DD[ind.tri] <- new.low.tri.vals
               DD <- t(DD)
               DD[ind.tri] <- new.low.tri.vals
          }

          attr(DD, "class") <- "crossdist"
          attr(DD, "method") <- "LB_Keogh1"

     } else {

          if (is.matrix(y))
               y <- lapply(seq_len(nrow(y)), function(i) y[i,])
          else if (is.numeric(y))
               y <- list(y)
          else if (is.list(y))
               y <- y
          else
               stop("Unsupported type for y")

          if (error.check)
               consistency_check(y, "tslist")

          if (window.size > length(y[[1]]))
               stop("Window size should not exceed length of the time series")

          # from 'caTools' package
          ## NOTE: the 'window.size' definition varies betwen 'dtw' and 'runmax/min'
          upper.env <- lapply(y, runmax, k=window.size*2+1, endrule="constant")
          lower.env <- lapply(y, runmin, k=window.size*2+1, endrule="constant")

          DD <- sapply(X=x, U=upper.env, L=lower.env,
                       FUN = function(x, ...) {
                            U <- list(...)$U
                            L <- list(...)$L

                            ## This will return one column of the distance matrix
                            D <- mapply(U, L, MoreArgs=list(x=x),
                                        FUN = function(u, l, x) {

                                             D <- rep(0, length(x))

                                             ind1 <- which(x > u)
                                             D[ind1] <- x[ind1] - u[ind1]
                                             ind2 <- which(x < l)
                                             D[ind2] <- l[ind2] - x[ind2]

                                             d <- switch(EXPR = norm,
                                                         L1 = sum(D),
                                                         L2 = sqrt(sum(D^2)))

                                             d

                                        })

                            D
                       })

          if (force.symmetry) {
               ind.tri <- lower.tri(DD)

               new.low.tri.vals <- t(DD)[ind.tri]
               indCorrect <- DD[ind.tri] > new.low.tri.vals
               new.low.tri.vals[indCorrect] <- DD[ind.tri][indCorrect]

               DD[ind.tri] <- new.low.tri.vals
               DD <- t(DD)
               DD[ind.tri] <- new.low.tri.vals
          }

          attr(DD, "class") <- "crossdist"
          attr(DD, "method") <- "LB_Keogh1"
     }

     t(DD)
}
