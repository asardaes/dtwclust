#' Wrapper for z-normalization
#'
#' Wrapper for function [base::scale()] that returns zeros instead of `NaN`. It also supports
#' matrices, data frames, and lists of time series.
#'
#' @export
#'
#' @param x Data to normalize. Either a vector, a matrix/data.frame where each row is to be
#'   normalized, or a list of vectors/matrices.
#' @param ... Further arguments to pass to [base::scale()].
#' @param multivariate Is `x` a multivariate time series? It will be detected automatically if a
#'   list is provided in `x`.
#' @param keep.attributes Should the mean and standard deviation returned by [base::scale()] be
#'   preserved?
#' @param error.check `r roxygen_error_check_param()`
#'
#' @details
#'
#' Multivariate series must have time spanning the rows and variables spanning the columns.
#'
#' @return Normalized data in the same format as provided.
#'
zscore <- function(x, ..., multivariate = FALSE, keep.attributes = FALSE, error.check = TRUE) {
    if (is.list(x) && !is.data.frame(x)) {
        x <- lapply(x, zscore, ...,
                    multivariate = is_multivariate(x),
                    keep.attributes = keep.attributes)
    }
    else if (!multivariate && (is.matrix(x) || is.data.frame(x))) {
        if (is.data.frame(x)) x <- base::as.matrix(x)
        if (error.check) check_consistency(x, "ts")
        dots <- list(...)
        center <- if (is.null(dots$center)) formals(base::scale)$center else dots$center
        scale <- if (is.null(dots$scale)) formals(base::scale)$scale else dots$scale
        x <- t(base::scale(t(x), center = center, scale = scale))
        x[is.nan(x)] <- 0
        if (!keep.attributes) { attr(x, "scaled:center") <- attr(x, "scaled:scale") <- NULL }
    }
    else {
        if (is.data.frame(x)) x <- base::as.matrix(x)
        if (error.check) check_consistency(x, "ts")
        dots <- list(...)
        center <- if (is.null(dots$center)) formals(base::scale)$center else dots$center
        scale <- if (is.null(dots$scale)) formals(base::scale)$scale else dots$scale
        x <- base::scale(x, center = center, scale = scale)
        x[is.nan(x)] <- 0
        if (!multivariate) dim(x) <- NULL
        if (!keep.attributes) { attr(x, "scaled:center") <- attr(x, "scaled:scale") <- NULL }
    }
    # return
    x
}
