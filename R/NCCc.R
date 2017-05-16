#' Cross-correlation with coefficient normalization
#'
#' This function uses FFT to compute the cross-correlation sequence between two series. They need
#' not be of equal length.
#'
#' @export
#'
#' @param x,y Univariate time series.
#'
#' @return The cross-correlation sequence with length `length(x) + length(y) - 1L`.
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In *Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data*, series
#' SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' @seealso
#'
#' [SBD()]
#'
NCCc <- function(x, y) {
    den <- lnorm(x, 2) * lnorm(y, 2)
    # Notice that the native 'convolve' function already uses FFT for the calculation
    if (den == 0) Inf else { stats::convolve(x, y, conj = TRUE, type = "open") / den }
}
