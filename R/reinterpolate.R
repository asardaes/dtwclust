#' Wrapper for simple linear reinterpolation
#'
#' This function is just a wrapper for the native function \code{\link[stats]{approx}} to do simple
#' linear reinterpolation. It also supports matrices, data frames, and lists of time series.
#'
#' @export
#'
#' @param x Data to reinterpolate. Either a vector, a matrix/data.frame where each row is to be
#'   reinterpolated, or a list of vectors/matrices.
#' @param new.length Desired length of the output series.
#' @param multivariate Is \code{x} a multivariate time series? It will be detected automatically if
#'   a list is provided in \code{x}.
#'
#' @details
#'
#' Multivariate series must have time spanning the rows and variables spanning the columns.
#'
#' @return Reinterpolated time series
#'
#' @examples
#'
#' data(uciCT)
#'
#' # list of univariate series
#' series <- reinterpolate(CharTraj, 205L)
#'
#' # list of multivariate series
#' series <- reinterpolate(CharTrajMV, 205L)
#'
#' # single multivariate series
#' series <- reinterpolate(CharTrajMV[[1L]], 205L, TRUE)
#'
reinterpolate <- function(x, new.length, multivariate = FALSE) {
    if (is.list(x) && !is.data.frame(x)) {
        x <- lapply(x, reinterpolate, new.length = new.length, multivariate = is_multivariate(x))

    } else if (!multivariate && (is.matrix(x) || is.data.frame(x))) {
        x <- t(apply(x, 1L, reinterpolate, new.length = new.length))

    } else {
        check_consistency(x, "ts")

        if (multivariate)
            x <- apply(x, 2L, reinterpolate, new.length = new.length)
        else
            x <- stats::approx(x, method = "linear", n = new.length)$y
    }

    x
}
