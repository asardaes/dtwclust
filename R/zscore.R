#' Wrapper for z-normalization
#'
#' Wrapper for function \code{\link[base]{scale}} that returns zeros instead of \code{NaN}. It also
#' supports a list of vectors and a matrix input.
#'
#' @param x Data to normalize. Either a vector, a matrix/data.frame where each row is to be normalized, or a list of
#' vectors.
#' @param ... Further arguments to pass to \code{\link[base]{scale}}.
#' @param na.rm Logical flag. Should \code{NA}s be removed? Ignored for matrix input.
#' @param multivariate Is \code{x} a multivariate time series? It will be detected automatically if a list is provided in
#' \code{x}.
#'
#' @return Normalized data in the same format as provided.
#'
#' @export
#'

zscore <- function(x, ..., na.rm = FALSE, multivariate = FALSE) {
     if (is.list(x)) {
          x <- lapply(x, zscore, na.rm = na.rm, multivariate = !is.null(dim(x[[1L]])), ...)

     } else if (!multivariate && (is.matrix(x) || is.data.frame(x))) {
          x <- t(apply(x, 1L, zscore, na.rm = na.rm, ...))

     } else {
          consistency_check(x, "ts")

          if (multivariate) {
               x <- apply(x, 2L, zscore, na.rm = na.rm, ...)

          } else {
               if (na.rm) x <- x[!is.na(x)]

               x <- scale(x, ...)
               x <- as.numeric(x) # scale returns columns
               x[is.nan(x)] <- 0
          }
     }

     x
}
