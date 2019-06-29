#' Shape average of several time series
#'
#' Time-series shape extraction based on optimal alignments as proposed by Paparrizos and Gravano
#' (2015) for the k-Shape clustering algorithm.
#'
#' @export
#' @importFrom RSpectra eigs_sym
#'
#' @param X A matrix or data frame where each row is a time series, or a list where each element is
#'   a time series. Multivariate series should be provided as a list of matrices where time spans
#'   the rows and the variables span the columns.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of `X`
#'   if `NULL`. For multivariate series, this should be a matrix with the same characteristics as
#'   the matrices in `X`. **It will be z-normalized**.
#' @param znorm Logical flag. Should z-scores be calculated for `X` before processing?
#' @param ... Further arguments for [zscore()].
#' @template error-check
#'
#' @details
#'
#' This works only if the series are *z-normalized*, since the output will also have this
#' normalization.
#'
#' The resulting centroid will have the same length as `centroid` if provided. Otherwise, there are
#' two possibilities: if all series from `X` have the same length, all of them will be used as-is,
#' and the output will have the same length as the series; if series have different lengths, a
#' series will be chosen at random and used as reference. The output series will then have the same
#' length as the chosen series.
#'
#' This centroid computation is cast as an optimization problem called maximization of Rayleigh
#' Quotient. It depends on the [SBD()] algorithm. See the cited article for more details.
#'
#' @return Centroid time series (z-normalized).
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In *Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data*, series
#' SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' @seealso
#'
#' [SBD()], [zscore()]
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
shape_extraction <- function(X, centroid = NULL, znorm = FALSE, ..., error.check = TRUE) {
    X <- tslist(X)
    if (error.check) {
        check_consistency(X, "vltslist")
        if (!is.null(centroid)) check_consistency(centroid, "ts")
    }

    # UTILS-utils.R
    if (is_multivariate(X)) {
        if (!is.null(centroid) && ncol(X[[1L]]) != NCOL(centroid))
            stop("Dimension inconsistency between the series in 'X' and the provided 'centroid'.")
        mv <- reshape_multivariate(X, centroid) # UTILS-utils.R
        new_c <- Map(mv$series, mv$cent, f = function(xx, cc, ...) {
            new_c <- shape_extraction(xx, cc, znorm = znorm, ..., error.check = FALSE)
        })
        return(do.call(cbind, new_c, TRUE))
    }

    Xz <- if (znorm) zscore(X, ..., error.check = FALSE) else X
    # make sure at least one series is not just a flat line at zero
    if (all(sapply(Xz, sum) == 0)) {
        if (is.null(centroid))
            return(rep(0, sample(lengths(Xz), 1L)))
        else
            return(centroid)
    }

    if (is.null(centroid)) {
        if (!different_lengths(Xz)) {
            A <- do.call(rbind, Xz, TRUE) # use all
        }
        else {
            centroid <- Xz[[sample(length(Xz), 1L)]] # random choice as reference
            A <- lapply(Xz, function(a) { SBD(centroid, a)$yshift })
            A <- do.call(rbind, A, TRUE)
        }
    }
    else {
        centroid <- zscore(centroid, ..., error.check = FALSE) # use given reference
        A <- lapply(Xz, function(a) { SBD(centroid, a)$yshift })
        A <- do.call(rbind, A, TRUE)
    }

    Y <- zscore(A, ..., error.check = FALSE)
    S <- if (is.matrix(Y)) t(Y) %*% Y else Y %*% t(Y)
    nc <- ncol(A)
    P <- diag(nc) - 1 / nc * matrix(1, nc, nc)
    M <- P %*% S %*% P
    ksc <- Re(RSpectra::eigs_sym(M, 1L)$vectors[ , 1L, drop = TRUE])
    # UTILS-utils.R
    d1 <- l2norm(A[1L, , drop = TRUE] - ksc)
    d2 <- l2norm(A[1L, , drop = TRUE] + ksc)
    if (d1 >= d2) ksc <- -ksc
    ksc <- zscore(ksc, ..., error.check = FALSE)
    ksc
}
