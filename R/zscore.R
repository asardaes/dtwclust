#' Wrapper for z-normalization
#'
#' Wrapper for function \code{\link[base]{scale}} that returns zeros instead of \code{NaN}. It also
#' supports a list of vectors and a matrix input.
#'
#' @param x Data to normalize. Either a vector, a matrix/data.frame where each row is to be normalized, or a list of
#' vectors.
#' @param ... Further arguments to pass to \code{\link[base]{scale}}.
#' @param multivariate Is \code{x} a multivariate time series? It will be detected automatically if a list is provided in
#' \code{x}.
#' @param na.rm Deprecated
#'
#' @return Normalized data in the same format as provided.
#'
#' @export
#'

zscore <- function(x, ..., na.rm, multivariate = FALSE) {
     if (!missing(na.rm))
          warning("The 'na.rm' has been deprecated.")

     if (is.list(x)) {
          x <- lapply(x, zscore, multivariate = !is.null(dim(x[[1L]])), ...)

     } else if (!multivariate && (is.matrix(x) || is.data.frame(x))) {
          x <- t(apply(x, 1L, zscore, ...))

     } else {
          consistency_check(x, "ts")

          dots <- list(...)
          center <- if(is.null(dots$center)) formals(scale)$center else dots$center
          scale <- if(is.null(dots$scale)) formals(scale)$scale else dots$scale

          x <- scale(x, center = center, scale = scale)
          x[is.nan(x)] <- 0

          if (multivariate)
               attr(x, "scaled:center") <- attr(x, "scaled:scale") <- NULL
          else
               x <- as.numeric(x) # remove dimension and other attributes
     }

     x
}
