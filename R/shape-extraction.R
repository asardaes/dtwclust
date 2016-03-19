#' Shape average of several time series
#'
#' Time-series shape extraction based on optimal alignments as proposed by Paparrizos and Gravano, 2015, for
#' the k-Shape clustering algorithm.
#'
#' This works only if the series are \emph{z-normalized}, since the output will also have this normalization.
#'
#' The resulting centroid will have the same length as \code{center} if provided. Otherwise, there are two
#' possibilities: if all series from \code{X} have the same length, all of them
#' will be used as-is, and the output will have the same length as the series; if series have different
#' lengths, a series will be chosen at random and used as reference. The output series will then have the
#' same length as the chosen series.
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
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In \emph{Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data},
#' series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
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
#' @param center Center to use as basis. \emph{It will be z-normalized}.
#' @param znorm Logical flag. Should z-scores be calculated for \code{X} before processing?
#'
#' @return Centroid time series.
#'
#' @export
#'

shape_extraction <- function(X, center = NULL, znorm = FALSE) {

     X <- consistency_check(X, "tsmat")

     consistency_check(X, "vltslist")

     if (znorm)
          Xz <- zscore(X)
     else
          Xz <- X

     ## make sure at least one series is not just a flat line at zero
     if (all(sapply(Xz, sum) == 0)) {
          if (is.null(center)) {
               return(rep(0, sample(lengths(Xz), 1)))

          } else {
               return(center)
          }
     }

     if (is.null(center)) {
          if (length(unique(lengths(Xz))) == 1L)
               A <- do.call(rbind, Xz) # use all
          else {
               center <- Xz[[sample(length(Xz), 1L)]] # random choice as reference

               A <- lapply(Xz, function(a) { SBD(center, a)$yshift })

               A <- do.call(rbind, A)
          }

     } else {
          center <- zscore(center) # use given reference

          A <- lapply(Xz, function(a) { SBD(center, a)$yshift })

          A <- do.call(rbind, A)
     }

     Y <- zscore(A)

     if (is.matrix(Y))
          S <- t(Y) %*% Y
     else
          S <- Y %*% t(Y)

     nc <- ncol(A)
     P <- diag(nc) - 1 / nc * matrix(1, nc, nc)
     M <- P %*% S %*% P

     ksc <- eigen(M)$vectors[ , 1L]

     d1 <- sqrt(crossprod(A[1L, ] - ksc))
     d2 <- sqrt(crossprod(A[1L, ] + ksc))

     if (d1 >= d2)
          ksc <- -ksc

     ksc <- zscore(ksc)

     ksc
}
