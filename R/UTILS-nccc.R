#' Cross-correlation with coefficient normalization
#'
#' This function uses the FFT to compute the cross-correlation sequence between two series. They
#' need *not* be of equal length.
#'
#' @export
#' @importFrom stats convolve
#'
#' @param x,y Univariate time series.
#' @param error.check `r roxygen_error_check_param()`
#'
#' @return The cross-correlation sequence with length `length(x) + length(y) - 1L`.
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In *Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data*, series
#' SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \doi{10.1145/2723372.2737793}.
#'
#' @seealso
#'
#' [SBD()]
#'
NCCc <- function(x, y, error.check = TRUE) {
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
        if (is_multivariate(list(x,y))) stop("NCCc does not support multivariate series.")
        x <- as.numeric(x)
        y <- as.numeric(y)
    }
    den <- l2norm(x) * l2norm(y) # UTILS-utils.R
    # Notice that the native 'convolve' function already uses FFT for the calculation
    if (den == 0) Inf else { stats::convolve(x, y, conj = TRUE, type = "open") / den }
}
