#' Shape average of several time series
#'
#' Time-series shape extraction based on optimal alignments as proposed by Papparizos and Gravano, 2015, for
#' the k-Shape clustering algorithm.
#'
#' This works only if the signals are \emph{z-normalized}, since the output will also have this normalization.
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
#' # Subset of interest, normalized
#' X <- t(sapply(CharTraj[1:5], zscore))
#'
#' # Obtain centroid series
#' C <- shape_extraction(X, znorm = FALSE)
#'
#' # Result
#' matplot(t(X), type = "l", col = 1:5)
#' points(C)
#'
#' @param X Numeric matrix where each row is a time series.
#' @param cz Center to use as basis. It should already be \emph{normalized}. Calculation uses all \code{X}
#' if \code{cz = NULL}.
#' @param znorm Boolean flag. Should z-scores be calculated for \code{X} before processing?
#'
#' @return Centroid time series.
#'
#' @export

shape_extraction <- function(X, cz = NULL, znorm = FALSE) {

     if (!is.matrix(X))
          stop("Unsupported type for X")

     if (znorm)
          Xz <- t(apply(X, 1, zscore))
     else
          Xz <- X

     if (!is.null(cz)) {

          cz <- as.numeric(cz)

          if (all(cz == 0)) {
               a <- Xz

          } else {
               a <- t(apply(Xz, 1, function(A) {
                    sbd <- SBD(cz, A)

                    sbd$yshift
               }))
          }

     } else {

          a <- Xz

     }

     Y <- t(apply(a, 1, zscore))

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
     attributes(ksc) <- NULL
     ksc <- as.numeric(ksc)

     ksc
}
