#' Wrapper for z-normalization
#'
#' Wrapper for function \code{\link[base]{scale}} that returns zeros instead of \code{NaN}.
#'
#' @param x Data to normalize.
#' @param ... Further arguments to pass to \code{\link[base]{scale}}.
#'
#' @return Normalized data.
#'
#' @export

zscore <- function(x, ...) {
     x <- scale(x, ...)
     x[is.nan(x)] <- 0
     dim(x) <- NULL # scale returns columns

     x
}
