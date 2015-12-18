#' Wrapper for simple linear reinterpolation
#'
#' This function is just a wrapper for the native function \code{\link[stats]{approx}} to do simple linear reinterpolation.
#'
#' @param ts A time series.
#' @param newLength Desired length of the output series.
#'
#' @return Reinterpolated time series
#'
#' @export
#'

reinterpolate <- function(ts, newLength) {

     newTS <- stats::approx(ts, method='linear', n=newLength)

     newTS$y
}
