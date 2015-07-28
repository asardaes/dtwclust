#' Keogh's lower bound
#'
#' This function calculates a lower bound (LB) on the Dynamic Time Warp (DTW) distance between two time
#' series. It uses a Sakoe-Chiba constraint.
#'
#' Because of the way the different functions being used here are implemented, there is a subtle but critical
#' mismatch in the way the window size is defined for DTW and the LB. The LB calculation expects an \emph{odd}
#' \code{window.size} that represents the whole window width to be used in the running max and min. The
#' outcome of said running functions are centered with respect to the window. The DTW calculation with
#' \code{\link[dtw]{dtw}} expects a window.size that represents the distance between the diagonal and one of
#' the edges of the window.
#'
#' Therefore, if, for example, the LB is calculated with a window of 21, the corresponding DTW distance should
#' be calculated with \code{21 \%/\% 2 = 10}.
#'
#' The internal functions take care of this discrepancy if needed, but you should be careful if you are
#' testing things manually.
#'
#' @param x A time series.
#' @param y A time series with the same length as \code{x}.
#' @param window.size Window size for envelop calculation. \strong{See details}.
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

lb_keogh <- function(x, y, window.size, norm = "L1") {

     norm <- match.arg(norm, c("L1", "L2"))

     consistency_check(x, "ts")
     consistency_check(y, "ts")

     if (length(x) != length(y)) {
          stop("The series must have the same length")
     }
     if (window.size%%2 != 1) {
          stop("For the Sakoe-Chiba band, the window must be symmetric and window.size must be odd")
     }
     if (window.size > length(x)) {
          stop("The width of the window should not exceed the length of the series")
     } else if (window.size <= 1) {
          stop("Window width must be larger than 1")
     }

     # from 'caTools' package
     upper.env <- runmax(y, window.size)
     lower.env <- runmin(y, window.size)

     D <- rep(0, length(x))

     # Using L1 in accordance with Lemire 2009
     ind1 <- which(x > upper.env)
     D[ind1] <- x[ind1] - upper.env[ind1]
     ind2 <- which(x < lower.env)
     D[ind2] <- lower.env[ind2] - x[ind2]

     d <- switch(EXPR = norm,
                 L1 = sum(D),
                 L2 = sqrt(sum(D^2)))

     return(list(d = d,
                 upper.env = upper.env,
                 lower.env = lower.env))
}

# ========================================================================================================
# Loop without using native 'proxy' (to avoid multiple calculations of the envelop)
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


     if (is.null(window.size)) {
          stop("Please provide the 'window.size' parameter")
     }
     if (window.size%%2 != 1) {
          stop("For the Sakoe-Chiba band, the window must be symmetric and window.size must be odd")
     }
     if (window.size <= 1) {
          stop("Window width must be larger than 1")
     }

     ## For looping convenience
     if (class(x) == "matrix")
          x <- lapply(seq_len(nrow(x)), function(i) x[i,])
     if (class(y) == "matrix")
          y <- lapply(seq_len(nrow(y)), function(i) y[i,])

     if (error.check)
          consistency_check(x, "tslist")

     if (window.size > length(x[[1]]))
          stop("Window size should not exceed length of the time series")

     if (is.null(y)) {
          # from 'caTools' package
          upper.env <- lapply(x, runmax, k=window.size, endrule="constant")
          lower.env <- lapply(x, runmin, k=window.size, endrule="constant")

          DD <- sapply(X=x, U=upper.env, L=lower.env,
                       FUN = function(x, ...) {
                            U <- list(...)$U
                            L <- list(...)$L

                            ## This will return one column of the distance matrix
                            D <- mapply(U, L, MoreArgs=list(x=x),
                                        FUN = function(u, l, x) {

                                             D <- rep(0, length(x))

                                             # Using L1 in accordance with Lemire 2009
                                             ind1 <- which(x > u)
                                             D[ind1] <- x[ind1] - u[ind1]
                                             ind2 <- which(x < l)
                                             D[ind2] <- l[ind2] - x[ind2]

                                             d <- switch(EXPR = norm,
                                                         L1 = sum(D),
                                                         L2 = sqrt(sum(D^2)))

                                             return(d)

                                        })

                            return(D)
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
          return(t(DD))


     } else {

          if (error.check)
               consistency_check(y, "tslist")

          if (window.size > length(y[[1]]))
               stop("Window size should not exceed length of the time series")

          # from 'caTools' package
          upper.env <- lapply(y, runmax, k=window.size, endrule="constant")
          lower.env <- lapply(y, runmin, k=window.size, endrule="constant")

          DD <- sapply(X=x, U=upper.env, L=lower.env,
                       FUN = function(x, ...) {
                            U <- list(...)$U
                            L <- list(...)$L

                            ## This will return one column of the distance matrix
                            D <- mapply(U, L, MoreArgs=list(x=x),
                                        FUN = function(u, l, x) {

                                             D <- rep(0, length(x))

                                             # Using L1 in accordance with Lemire 2009
                                             ind1 <- which(x > u)
                                             D[ind1] <- x[ind1] - u[ind1]
                                             ind2 <- which(x < l)
                                             D[ind2] <- l[ind2] - x[ind2]

                                             d <- switch(EXPR = norm,
                                                         L1 = sum(D),
                                                         L2 = sqrt(sum(D^2)))

                                             return(d)

                                        })

                            return(D)
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
          return(t(DD))
     }
}
