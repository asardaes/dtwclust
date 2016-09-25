#' Shape average of several time series
#'
#' Time-series shape extraction based on optimal alignments as proposed by Paparrizos and Gravano, 2015, for
#' the k-Shape clustering algorithm.
#'
#' This works only if the series are \emph{z-normalized}, since the output will also have this normalization.
#'
#' The resulting centroid will have the same length as \code{centroid} if provided. Otherwise, there are two
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
#' @param X A data matrix where each row is a time series, or a list where each element is a time series.
#' Multivariate series should be provided as a list of matrices where time spans the rows and the variables
#' span the columns.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of \code{X} if
#' \code{NULL}. For multivariate series, this should be a matrix with the same characteristics as the
#' matrices in \code{X}. \emph{It will be z-normalized}.
#' @param center Deprecated, please use \code{centroid} instead.
#' @param znorm Logical flag. Should z-scores be calculated for \code{X} before processing?
#'
#' @return Centroid time series (z-normalized).
#'
#' @export
#'
shape_extraction <- function(X, centroid = NULL, center = NULL, znorm = FALSE) {
     if (!missing(center)) {
          warning("The 'center' argument has been deprecated and will be removed in the next version. ",
                  "Please use 'centroid' instead.")

          if (is.null(centroid)) centroid <- center
     }

     X <- consistency_check(X, "tsmat")

     consistency_check(X, "vltslist")

     ## utils.R
     if (check_multivariate(X)) {
          ## multivariate
          mv <- reshape_multviariate(X, centroid) # utils.R

          new_c <- mapply(mv$series, mv$cent, SIMPLIFY = FALSE,
                          FUN = function(xx, cc) {
                               new_c <- shape_extraction(xx, cc, znorm = znorm)
                          })

          return(do.call(cbind, new_c))
     }

     if (znorm)
          Xz <- zscore(X)
     else
          Xz <- X

     ## make sure at least one series is not just a flat line at zero
     if (all(sapply(Xz, sum) == 0)) {
          if (is.null(centroid)) {
               return(rep(0, sample(lengths(Xz), 1L)))

          } else {
               return(centroid)
          }
     }

     if (is.null(centroid)) {
          if (!check_lengths(Xz)) {
               A <- do.call(rbind, Xz) # use all

          } else {
               centroid <- Xz[[sample(length(Xz), 1L)]] # random choice as reference

               A <- lapply(Xz, function(a) { SBD(centroid, a)$yshift })

               A <- do.call(rbind, A)
          }

     } else {
          consistency_check(centroid, "ts")

          centroid <- zscore(centroid) # use given reference

          A <- lapply(Xz, function(a) { SBD(centroid, a)$yshift })

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

     ksc <- eigen(M)$vectors[ , 1L, drop = TRUE]

     d1 <- lnorm(A[1L, , drop = TRUE] - ksc, 2)
     d2 <- lnorm(A[1L, , drop = TRUE] + ksc, 2)

     if (d1 >= d2)
          ksc <- -ksc

     ksc <- zscore(ksc)

     ksc
}
