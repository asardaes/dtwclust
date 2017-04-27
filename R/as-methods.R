# ==================================================================================================
# Coercion methods for cross/pair-dist
# ==================================================================================================

#' as.matrix
#'
#' \pkg{proxy} exported a non-generic `as.matrix` function. This is to re-export the base version
#' and add some coercion methods for `pairdist` and `crossdist`.
#'
#' @name as.matrix
#' @export
#'
#' @param x,... See [base::as.matrix()].
#'
#' @seealso [base::as.matrix()]
#'
setGeneric("as.matrix", package = "base")

#' @method as.matrix crossdist
#' @export
#'
as.matrix.crossdist <- function(x, ...) {
    x <- cbind(x)
    class(x) <- "matrix"
    x
}

#' @method as.matrix pairdist
#' @export
#'
as.matrix.pairdist <- function(x, ...) {
    x <- cbind(x)
    class(x) <- "matrix"
    x
}

#' @method as.data.frame crossdist
#' @export
#'
as.data.frame.crossdist <- function(x, ...) {
    as.data.frame(as.matrix(x, ...), ...)
}

#' @method as.data.frame pairdist
#' @export
#'
as.data.frame.pairdist <- function(x, ...) {
    as.data.frame(as.matrix(x, ...), ...)
}
