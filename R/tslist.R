#' Coerce matrices or data frames to a list of time series
#'
#' Change a matrix or data frame to a list of univariate time series
#'
#' @export
#'
#' @param series A matrix or data frame where each row is a time series.
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
#' The function assumes that indexing the matrix/data.frame with `series[i, ]` will return something
#' that can be coerced to a numeric vector.
#'
#' No consistency checks are performed by this function.
#'
tslist <- function(series) {
    if (is.matrix(series) || is.data.frame(series)) {
        rnms <- rownames(series)
        series <- lapply(seq_len(nrow(series)), function(i) { as.numeric(series[i, ]) })
        if (!is.null(rnms)) setnames_inplace(series, rnms)

    } else if (is.numeric(series))
        series <- list(series)
    else if (!is.list(series))
        stop("Unsupported data type.")

    ## return
    series
}
