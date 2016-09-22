#' Basic DTW distance
#'
#' This is a custom implementation of the DTW algorithm without all the functionality included in
#' \code{\link[dtw]{dtw}}. Because of that, it should be slightly faster, while still supporting the most
#' common options.
#'
#' If \code{backtrack} is \code{TRUE}, the mapping of indices between series is returned in a list. The
#' backtracking algorithm is faster than the one in \code{\link[dtw]{dtw}}, but the results differ when more
#' than one optimal warping path exists.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10(2) + 1 = 21} observations falling within
#' the window.
#'
#' @param x,y Time series. Multivariate series must have time spanning the rows and variables spanning the
#' columns.
#' @param window.size Size for slanted band window. \code{NULL} means no constraint.
#' @param norm Norm for the DTW calculation, "L1" for Manhattan or "L2" for Euclidean.
#' @param step.pattern Step pattern for DTW. Only \code{symmetric1} or \code{symmetric2} supported here. See
#' \code{\link[dtw]{stepPattern}}.
#' @param backtrack Also compute the warping path between series? See details.
#' @param normalize Should the distance be normalized? Only supported for \code{symmetric2}.
#' @param ... Currently ignored.
#' @param gcm Optionally, a matrix with \code{NROW(x)+1} rows and \code{NROW(y)+1} columns to use for
#' the global cost matrix calculations. Used internally for memory optimization.
#' @param dm Optionally, a matrix with \code{NROW(x)+1} rows and \code{NROW(y)+1} columns to use for
#' the direction matrix and backtracking calculations. Used internally for memory optimization.
#'
#' @return The DTW distance. For \code{backtrack} \code{=} \code{TRUE}, a list with: \itemize{
#'   \item \code{distance}: The DTW distance.
#'   \item \code{index1}: \code{x} indices for the matched elements in the warping path.
#'   \item \code{index2}: \code{y} indices for the matched elements in the warping path.
#' }
#'
#' @export
#'
dtw_basic <- function(x, y, window.size = NULL, norm = "L1",
                      step.pattern = get("symmetric2"), backtrack = FALSE,
                      normalize = FALSE, ..., gcm = NULL, dm = NULL) {
     consistency_check(x, "ts")
     consistency_check(y, "ts")

     backtrack <- as.logical(backtrack)

     if (NCOL(x) != NCOL(y))
          stop("Multivariate series must have the same number of variables.")

     if (is.null(window.size))
          window.size <- -1L
     else
          window.size <- consistency_check(window.size, "window")

     norm <- match.arg(norm, c("L1", "L2"))
     norm <- switch(norm, "L1" = 1, "L2" = 2)

     if (identical(step.pattern, get("symmetric1")))
          step.pattern <- 1
     else if (identical(step.pattern, get("symmetric2")))
          step.pattern <- 2
     else
          stop("step.pattern must be either symmetric1 or symmetric2")

     if (is.null(gcm))
          gcm <- matrix(-1, NROW(x) + 1L, NROW(y) + 1L)
     else if (nrow(gcm) < NROW(x) + 1L || ncol(gcm) < NROW(y) + 1L)
          stop("dtw_basic: Dimension inconsistency in 'gcm'")

     if (backtrack) {
          if (is.null(dm))
               dm <- matrix(-1L, NROW(x) + 1L, NROW(y) + 1L)
          else {
               storage.mode(dm) <- "integer"

               if (nrow(dm) < NROW(x) + 1L || ncol(dm) < NROW(y) + 1L)
                    stop("dtw_basic: Dimension inconsistency in 'dm'")
          }
     }

     d <- .Call("dtw_basic", x, y, window.size,
                NROW(x), NROW(y), NCOL(x),
                norm, step.pattern, backtrack,
                gcm, dm,
                PACKAGE = "dtwclust")

     if (normalize && step.pattern == 2) {
          if (backtrack)
               d$distance <- d$distance / (NROW(x) + NROW(y))
          else
               d <- d / (NROW(x) + NROW(y))

     } else if (normalize && step.pattern != 2)
          warning("Unable to normalize with the chosen 'step.pattern'.")

     if (backtrack) {
          d$index1 <- d$index1[d$path:1L]
          d$index2 <- d$index2[d$path:1L]
          d$path <- NULL
     }

     d
}

dtw_basic_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1",
                            step.pattern = get("symmetric2"), backtrack = FALSE,
                            normalize = FALSE, ..., gcm = NULL, dm = NULL, pairwise = FALSE) {
     x <- consistency_check(x, "tsmat")
     consistency_check(x, "vltslist")

     if (is.null(y)) {
          y <- x

     } else {
          y <- consistency_check(y, "tsmat")
          consistency_check(y, "vltslist")
     }

     ## Register doSEQ if necessary
     check_parallel()

     X <- split_parallel(x)

     if (pairwise)
          Y <- split_parallel(y)
     else
          Y <- y

     ## Calculate distance matrix
     if (pairwise) {
          D <- foreach(x = X, y = Y,
                       .combine = c,
                       .multicombine = TRUE,
                       .packages = "dtwclust") %dopar% {
                            L1 <- max(lengths(x))
                            L2 <- max(lengths(y))

                            gcm <- matrix(0, L1 + 1, L2 + 1)
                            dm <- matrix(0L, L1 + 1, L2 + 1)

                            mapply(x, y, FUN = function(x, y) {
                                 dtw_basic(x, y, window.size = window.size,
                                           norm = norm, step.pattern = step.pattern,
                                           backtrack = FALSE, normalize = normalize,
                                           gcm = gcm, dm = dm)
                            })
                       }

          names(D) <- NULL
          attr(D, "class") <- "pairdist"

     } else {
          D <- foreach(x = X,
                       .combine = rbind,
                       .multicombine = TRUE,
                       .packages = "dtwclust") %dopar% {
                            L1 <- max(lengths(x))
                            L2 <- max(lengths(y))

                            gcm <- matrix(0, L1 + 1, L2 + 1)
                            dm <- matrix(0L, L1 + 1, L2 + 1)

                            ret <- lapply(x, y = Y, FUN = function(x, y) {
                                 sapply(y, x = x, FUN = function(y, x) {
                                      dtw_basic(x, y, window.size = window.size,
                                                norm = norm, step.pattern = step.pattern,
                                                backtrack = FALSE, normalize = normalize,
                                                gcm = gcm, dm = dm)
                                 })
                            })

                            do.call(rbind, ret)
                       }

          attr(D, "class") <- "crossdist"
     }

     attr(D, "method") <- "DTW_BASIC"

     D
}
