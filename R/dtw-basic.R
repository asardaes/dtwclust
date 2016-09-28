#' Basic DTW distance
#'
#' This is a custom implementation of the DTW algorithm without all the functionality included in
#' \code{\link[dtw]{dtw}}. Because of that, it should be slightly faster, while still supporting the most
#' common options.
#'
#' If \code{backtrack} is \code{TRUE}, the mapping of indices between series is returned in a list.
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
#' the global cost matrix calculations. Used internally for memory optimization. If provided, it \strong{will}
#' be modified \emph{in place} by \code{C} code, except in the parallel version in \code{proxy::}\code{\link[proxy]{dist}}
#' which ignores it for thread-safe reasons.
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
                      step.pattern = symmetric2, backtrack = FALSE,
                      normalize = FALSE, ..., gcm = NULL) {
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

     if (identical(step.pattern, symmetric1))
          step.pattern <- 1
     else if (identical(step.pattern, symmetric2))
          step.pattern <- 2
     else
          stop("step.pattern must be either symmetric1 or symmetric2")

     if (normalize && step.pattern == 1)
          stop("Unable to normalize with chosen step pattern.")

     if (is.null(gcm))
          gcm <- matrix(-1, NROW(x) + 1L, NROW(y) + 1L)
     else if (!is.matrix(gcm) || nrow(gcm) < (NROW(x) + 1L) || ncol(gcm) < (NROW(y) + 1L))
          stop("dtw_basic: Dimension inconsistency in 'gcm'")
     else if (storage.mode(gcm) != "double")
          stop("If provided, 'gcm' must have storage mode double.")

     d <- .Call("dtw_basic", x, y, window.size,
                NROW(x), NROW(y), NCOL(x),
                norm, step.pattern, backtrack,
                gcm, PACKAGE = "dtwclust")

     if (normalize && step.pattern == 2) {
          if (backtrack)
               d$distance <- d$distance / (NROW(x) + NROW(y))
          else
               d <- d / (NROW(x) + NROW(y))
     }

     if (backtrack) {
          d$index1 <- d$index1[d$path:1L]
          d$index2 <- d$index2[d$path:1L]
          d$path <- NULL
     }

     d
}

dtw_basic_proxy <- function(x, y = NULL, ..., gcm = NULL, pairwise = FALSE, symmetric = FALSE) {
     x <- consistency_check(x, "tsmat")
     consistency_check(x, "vltslist")

     dots <- list(...)

     if (is.null(y)) {
          y <- x
          symmetric <- is.null(dots$window.size) || !check_lengths(x)

     } else {
          y <- consistency_check(y, "tsmat")
          consistency_check(y, "vltslist")

          if (symmetric && length(x) != length(y))
               stop("dtw_basic: Dimension inconsistency when calculating cross-distance matrix within proxy.")
     }

     retclass <- "crossdist"

     if (symmetric && !pairwise) {
          check_parallel()

          pairs <- call_pairs(length(x), lower = FALSE)
          pairs <- split_parallel(pairs, 1L)

          dots$pairwise <- TRUE

          d <- foreach(pairs = pairs,
                       .combine = c,
                       .multicombine = TRUE,
                       .packages = "dtwclust",
                       .export = "enlist") %dopar% {
                            do.call(proxy::dist,
                                    enlist(x = x[pairs[ , 1L]],
                                           y = x[pairs[ , 2L]],
                                           method = "dtw_basic",
                                           dots = dots))
                       }

          rm("pairs")

          D <- matrix(0, nrow = length(x), ncol = length(x))
          D[upper.tri(D)] <- d
          D <- t(D)
          D[upper.tri(D)] <- d

          class(D) <- retclass
          attr(D, "dimnames") <- list(names(x), names(x))
          attr(D, "method") <- "DTW_BASIC"

          return(D)
     }

     X <- split_parallel(x)

     ## Register doSEQ if necessary
     if (check_parallel())
          GCM <- lapply(1L:length(X), function(dummy) NULL)
     else
          GCM <- list(gcm)

     ## Calculate distance matrix
     if (pairwise) {
          Y <- split_parallel(y)

          if (length(lengths(X)) != length(lengths(Y)) || any(lengths(X) != lengths(Y)))
               stop("Pairwise distances require the same amount of series in 'x' and 'y'")

          D <- foreach(x = X, y = Y, gcm = GCM,
                       .combine = c,
                       .multicombine = TRUE,
                       .packages = "dtwclust",
                       .export = "enlist") %dopar% {
                            L1 <- max(lengths(x))
                            L2 <- max(lengths(y))

                            if (is.null(gcm))
                                 dots$gcm <- matrix(0, L1 + 1L, L2 + 1L)
                            else if (!is.matrix(gcm) || nrow(gcm) < (L1 + 1L) || ncol(gcm) < (L2 + 1L))
                                 stop("dtw_basic: Dimension inconsistency in 'gcm'")
                            else if (storage.mode(gcm) != "double")
                                 stop("If provided, 'gcm' must have storage mode double.")

                            mapply(x, y, FUN = function(x, y) {
                                 do.call("dtw_basic",
                                         enlist(x = x,
                                                y = y,
                                                dots = dots))
                            })
                       }

          names(D) <- NULL
          retclass <- "pairdist"

     } else {
          D <- foreach(x = X, gcm = GCM,
                       .combine = rbind,
                       .multicombine = TRUE,
                       .packages = "dtwclust",
                       .export = "enlist") %dopar% {
                            L1 <- max(lengths(x))
                            L2 <- max(lengths(y))

                            if (is.null(gcm))
                                 dots$gcm <- matrix(0, L1 + 1, L2 + 1)
                            else if (!is.matrix(gcm) || nrow(gcm) < (L1 + 1L) || ncol(gcm) < (L2 + 1L))
                                 stop("dtw_basic: Dimension inconsistency in 'gcm'")
                            else if (storage.mode(gcm) != "double")
                                 stop("If provided, 'gcm' must have storage mode double.")

                            ret <- lapply(x, y = y, FUN = function(x, y) {
                                 sapply(y, x = x, FUN = function(y, x) {
                                      do.call("dtw_basic",
                                              enlist(x = x,
                                                     y = y,
                                                     dots = dots))
                                 })
                            })

                            do.call(rbind, ret)
                       }
     }

     class(D) <- retclass
     attr(D, "method") <- "DTW_BASIC"

     D
}
