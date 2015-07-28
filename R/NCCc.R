#' Cross-correlation with coefficient normalization
#'
#' This function uses FFT to compute the cross-correlation sequence between two series. They need not be of
#' equal length.
#'
#' @param x A series.
#' @param y Another series.
#'
#' @return The cross-correlation sequence with length \code{length(x) + length(y) - 1}.
#'
#' @export

NCCc <- function(x,y) {

     r <- convolve(x, y, type = "open")

     den <- sqrt(crossprod(x)) * sqrt(crossprod(y))

     CCseq <- r / den

     if (den == 0)
          return(Inf)
     else
          return(Re(CCseq))
}
