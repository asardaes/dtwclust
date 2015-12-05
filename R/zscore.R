#' Wrapper for z-normalization
#'
#' Wrapper for function \code{\link[base]{scale}} that returns zeros instead of \code{NaN}. It also
#' supports a list of vectors.
#'
#' @param x Data to normalize. Either a vector or a list of vectors.
#' @param ... Further arguments to pass to \code{\link[base]{scale}}.
#'
#' @return Normalized data.
#'
#' @export

zscore <- function(x, ...) {

     if (is.list(x)) {
          dots <- list(...)

          x <- lapply(x, function(xx) {
               xx <- do.call("scale", c(list(x=xx), dots))
               dim(xx) <- NULL # scale returns columns
               xx[is.nan(xx)] <- 0

               xx
          })

     } else {
          x <- scale(x, ...)
          dim(x) <- NULL # scale returns columns
          x[is.nan(x)] <- 0
     }

     x
}
