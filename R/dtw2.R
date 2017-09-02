#' DTW distance with L2 norm
#'
#' Wrapper for the [dtw::dtw()] function using L2 norm for both the local cost matrix (LCM) creation
#' as well as the final cost aggregation step.
#'
#' @export
#'
#' @param x,y A time series. A multivariate series should have time spanning the rows and variables
#'   spanning the columns.
#' @param ... Further arguments for [dtw::dtw()].
#'
#' @details
#'
#' The L-norms are used in two different steps by the DTW algorithm. First when creating the LCM,
#' where the element \eqn{(i,j)} of the matrix is computed as the L-norm of \eqn{x^v_i - y^v_j} for
#' all variables \eqn{v}. Note that this means that, in case of multivariate series, they must have
#' the same number of variables, and that univariate series will produce the same LCM regardless of
#' the L-norm used. After the warping path is found by DTW, the final distance is calculated as the
#' L-norm of all \eqn{(i,j)} elements of the LCM that fall on the warping path.
#'
#' The [dtw::dtw()] function allows changing the norm by means of its `dist.method` parameter, but
#' it only uses it when creating the LCM, and not when calculating the final aggregated cost, i.e.
#' the DTW distance.
#'
#' This wrapper simply returns the appropriate DTW distance using L2 norm (Euclidean norm). A
#' [proxy::dist()] version is also registered.
#'
#' @template window
#'
#' @return An object of class `dtw`.
#'
dtw2 <- function(x, y, ...) {
    lcm <- proxy::dist(x, y, method = "L1")
    d <- dtw::dtw(x = lcm^2, y = NULL, ...)
    d$distance <- sqrt(d$distance)

    if (!is.na(d$normalizedDistance)) {
        normalization <- switch(attr(d$stepPattern, "norm"),
                                "N" = nrow(lcm),
                                "M" = d$jmin,
                                "N+M" = nrow(lcm) + d$jmin,
                                stop("Unknown normalization factor for DTW."))

        d$normalizedDistance <- d$distance / normalization
    }

    ## return
    d
}

dtw2.proxy <- function(x, y, ...) {
    lcm <- proxy::dist(x, y, method = "L1")
    sqrt(dtw::dtw(x = lcm^2, y = NULL, distance.only = TRUE, ...)$distance)
}
