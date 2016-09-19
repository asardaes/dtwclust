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
                      step.pattern = get("symmetric2"), backtrack = FALSE, normalize = FALSE) {
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

     d <- .Call("dtw_basic", x, y, window.size,
                NROW(x), NROW(y), NCOL(x),
                norm, step.pattern, backtrack,
                PACKAGE = "dtwclust")

     if (normalize && step.pattern == 2)
          d <- d / (NROW(x) + NROW(y))

     if (backtrack) {
          path <- attr(d, "path")
          index1 <- attr(d, "index1")
          index1 <- index1[path:length(index1)]
          index2 <- attr(d, "index2")
          index2 <- index2[path:length(index2)]
          attributes(d) <- NULL

          d <- list(distance = d, index1 = index1, index2 = index2)
     }

     d
}

dtw_basic_proxy <- function(x, y, window.size = NULL, norm = "L1",
                            step.pattern = get("symmetric2"), backtrack = FALSE, normalize = FALSE) {
     dtw_basic(x, y, window.size = window.size,
               norm = norm, step.pattern = step.pattern,
               backtrack = FALSE, normalize = normalize)
}

dtw_dba <- function(x, y, window.size = NULL, norm = "L1",
                    step.pattern = get("symmetric2"), backtrack = TRUE,
                    normalize = FALSE, lcm_gcm = NULL, ...) {
     if (is.null(lcm_gcm))
          return(dtw_basic(x, y, window.size = window.size, norm = norm,
                           step.pattern = step.pattern, backtrack = TRUE,
                           normalize = normalize))
     else if (any(dim(lcm_gcm) != c(length(x) + 1L, length(y) + 1L)))
          stop("Dimension inconsistency when calculating DTW within DBA.")

     consistency_check(x, "ts")
     consistency_check(y, "ts")

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

     d <- .Call("dtw_dba", x, y, window.size,
                NROW(x), NROW(y), NCOL(x),
                norm, step.pattern, lcm_gcm,
                PACKAGE = "dtwclust")

     if (normalize && step.pattern == 2)
          d <- d / (NROW(x) + NROW(y))

     path <- attr(d, "path")
     index1 <- attr(d, "index1")
     index1 <- index1[path:length(index1)]
     index2 <- attr(d, "index2")
     index2 <- index2[path:length(index2)]
     attributes(d) <- NULL

     d <- list(distance = d, index1 = index1, index2 = index2)

     d
}
