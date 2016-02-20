#' Wrapper for z-normalization
#'
#' Wrapper for function \code{\link[base]{scale}} that returns zeros instead of \code{NaN}. It also
#' supports a list of vectors and a matrix input.
#'
#' @param x Data to normalize. Either a vector, a matrix where each row is to be normalized, or a list of
#' vectors.
#' @param ... Further arguments to pass to \code{\link[base]{scale}}.
#' @param na.rm Logical flag. Should \code{NA}s be removed? Ignored for matrix input.
#'
#' @return Normalized data in the same format as provided.
#'
#' @export
#'

zscore <- function(x, ..., na.rm = FALSE) {

     dots <- list(...)

     if (is.list(x)) {
          if (na.rm)
               x <- lapply(x, function(xx) { xx[!is.na(xx)] })

          consistency_check(x, "vltslist")

          x <- lapply(x, function(xx) {
               xx <- do.call("scale", c(list(x=xx), dots))
               xx <- as.numeric(xx) # scale returns columns
               xx[is.nan(xx)] <- 0

               xx
          })

     } else if (is.matrix(x)) {
          x <- t(apply(x, 1, function(xx) {
               xx <- do.call("scale", c(list(x=xx), dots))
               xx <- as.numeric(xx) # scale returns columns
               xx[is.nan(xx)] <- 0

               xx
          }))

     } else {
          if (na.rm)
               x <- x[!is.na(x)]

          consistency_check(x, "ts")

          x <- scale(x, ...)
          x <- as.numeric(x) # scale returns columns
          x[is.nan(x)] <- 0
     }

     x
}
