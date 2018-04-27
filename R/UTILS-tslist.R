#' Coerce matrices or data frames to a list of time series
#'
#' Change a matrix or data frame to a list of univariate time series
#'
#' @export
#'
#' @param series A matrix or data frame where each row is a time series.
#' @param simplify Coerce all series in the resulting list to either matrix (multivariate) or
#'   numeric (univariate).
#'
#' @details
#'
#' Almost all functions in \pkg{dtwclust} work internally with lists of time series. If you want to
#' avoid constant coercion, create a list of time series once by calling this function.
#'
#' For matrices and data frames, each **row** is considered as one time series. A list input is
#' simply passed through.
#'
#' @return
#'
#' A list of time series.
#'
#' @note
#'
#' The function assumes that matrix-like objects can be first coerced via [base::as.matrix()], so
#' that the result can be indexed with `series[i, ]`.
#'
#' No consistency checks are performed by this function.
#'
tslist <- function(series, simplify = FALSE) {
    if (is.matrix(series) || is.data.frame(series)) {
        rnms <- rownames(series)
        mat <- unname(base::as.matrix(series))
        series <- vector("list", nrow(mat))
        if (!is.null(rnms)) setnames_inplace(series, rnms)
        for (i in 1L:nrow(mat)) {
            series[[i]] <- mat[i,]
        }
    }
    else if (is.numeric(series))
        series <- list(series)
    else if (!is.list(series))
        stop("Unsupported data type.")
    # coerce to simple types that are known to work
    if (simplify) {
        if (is_multivariate(series))
            series <- lapply(series, base::as.matrix)
        else
            series <- lapply(series, base::as.numeric)
    }
    # return
    series
}
