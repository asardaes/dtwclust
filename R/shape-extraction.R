#' Shape average of several time series
#'
#' Time-series shape extraction based on optimal alignments as proposed by Papparizos and Gravano, 2015, for
#' the k-Shape clustering algorithm.
#'
#' This works only if the series are \emph{z-normalized}, since the output will also have this normalization.
#' The resulting centroid will have the same length as \code{cz} if provided, or the same length as the time
#' series in \code{X} if \code{cz = NULL}. Therefore, in the latter case, all series of \code{X} must have
#' equal lengths.
#'
#' This centroid computation is casted as an optimization problem called maximization of Rayleigh Quotient.
#' See the cited article for more details.
#'
#' @seealso
#'
#' \code{\link{SBD}}, \code{\link{zscore}}
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.'' In \emph{Proceedings of the 2015
#' ACM SIGMOD International Conference on Management of Data}, series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{
#' http://doi.org/10.1145/2723372.2737793}.
#'
#' @examples
#'
#' # Sample data
#' data(uciCT)
#'
#' # Normalize desired subset
#' X <- zscore(CharTraj[1:5])
#'
#' # Obtain centroid series
#' C <- shape_extraction(X)
#'
#' # Result
#' matplot(do.call(cbind, X),
#'         type = "l", col = 1:5)
#' points(C)
#'
#' @param X Numeric matrix where each row is a time series, or a list of time series.
#' @param cz Center to use as basis. \emph{It will be z-normalized}. Function uses all \code{X}
#' if \code{cz = NULL}. Must be provided and different than a zero-line if time series in
#' \code{X} have different lengths.
#' @param znorm Logical flag. Should z-scores be calculated for \code{X} before processing?
#'
#' @return Centroid time series.
#'
#' @export

shape_extraction <- function(X, cz = NULL, znorm = FALSE) {

     X <- consistency_check(X, "tsmat")

     if (is.null(cz))
          consistency_check(X, "tslist")
     else
          consistency_check(X, "vltslist")

     if (znorm)
          Xz <- zscore(X)
     else
          Xz <- X

     if (!is.null(cz)) {
          cz <- zscore(cz)

          if (all(cz == 0)) {
               a <- do.call(rbind, Xz)

          } else {
               a <- lapply(Xz, function(A) {
                    SBD(cz, A)$yshift
               })

               a <- do.call(rbind, a)
          }

     } else {
          a <- do.call(rbind, Xz)
     }

     Y <- zscore(a)

     if (is.matrix(Y))
          S <- t(Y) %*% Y
     else
          S <- Y %*% t(Y)

     nc <- ncol(a)
     P <- diag(nc) - 1 / nc * matrix(1, nc, nc)
     M <- P %*% S %*% P

     ksc <- eigen(M)$vectors[,1]

     d1 <- sqrt(crossprod(a[1,] - ksc))
     d2 <- sqrt(crossprod(a[1,] + ksc))

     if (d1 >= d2)
          ksc <- -ksc

     ksc <- zscore(ksc)

     ksc
}
