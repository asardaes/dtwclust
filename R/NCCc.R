#' Cross-correlation with coefficient normalization
#'
#' This function uses FFT to compute the cross-correlation sequence between two series. They need not be of
#' equal length.
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.'' In \emph{Proceedings of the 2015
#' ACM SIGMOD International Conference on Management of Data}, series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{
#' http://doi.org/10.1145/2723372.2737793}.
#'
#' @seealso
#'
#' \code{\link{SBD}}
#'
#' @param x A time series.
#' @param y Another time series.
#'
#' @return The cross-correlation sequence with length \code{length(x) + length(y) - 1}.
#'
#' @export
#' @importFrom stats convolve

NCCc <- function(x, y) {

     # Notice that the native 'convolve' function already uses FFT for the calculation
     r <- stats::convolve(x, y, conj = TRUE, type = "open")

     den <- sqrt(crossprod(x)) * sqrt(crossprod(y))

     CCseq <- r / den

     if (den == 0)
          return(Inf)
     else
          return(CCseq)
}
