#' Wrapper for simple linear reinterpolation
#'
#' This function is just a wrapper for the native function \code{\link[stats]{approx}} to do simple linear
#' reinterpolation.
#'
#' @param x Data to reinterpolate. Either a vector, a matrix/data.frame where each row is to be reinterpolated, or a list of
#' vectors.
#' @param new.length Desired length of the output series.
#' @param multivariate Is \code{x} a multivariate time series? It will be detected automatically if a list is provided in
#' \code{x}.
#' @param newLength Deprecated.
#'
#' @return Reinterpolated time series
#'
#' @export
#'

reinterpolate <- function(x, new.length, multivariate = FALSE, newLength) {
     if (!missing(newLength)) {
          warning("The 'newLength' argument has been deprecated, use 'new.length' instead.")

          if (missing(new.length)) new.length <- newLength
     }

     if (is.list(x)) {
          x <- lapply(x, reinterpolate, new.length = new.length, multivariate = !is.null(dim(x[[1L]])))

     } else if (!multivariate && (is.matrix(x) || is.data.frame(x))) {
          x <- t(apply(x, 1L, reinterpolate, new.length = new.length))

     } else {
          consistency_check(x, "ts")

          if (multivariate)
               x <- apply(x, 2L, reinterpolate, new.length = new.length)
          else
               x <- stats::approx(x, method = "linear", n = new.length)$y
     }

     x
}
