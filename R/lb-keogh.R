#' Keogh's DTW lower bound
#'
#' This function calculates a lower bound (LB) on the Dynamic Time Warp (DTW) distance between two
#' time series. It uses a Sakoe-Chiba constraint.
#'
#' @export
#'
#' @inheritParams lb_improved
#'
#' @details
#'
#' The reference time series should go in `x`, whereas the query time series should go in `y`.
#'
#' @template window
#'
#' @return A list with:
#'
#'   - `d`: The lower bound of the DTW distance.
#'   - `upper.env`: The time series of `y`'s upper envelop.
#'   - `lower.env`: The time series of `y`'s lower envelop.
#'
#' @note
#'
#' The lower bound is defined for time series of equal length only and is **not** symmetric.
#'
#' If you wish to calculate the lower bound between several time series, it would be better to use
#' the version registered with the `proxy` package, since it includes some small optimizations. The
#' convention mentioned above for references and queries still holds. See the examples.
#'
#' The proxy version of `force.symmetry` should only be used when only `x` is provided or both `x`
#' and `y` are identical. It compares the lower and upper triangular of the resulting distance
#' matrix and forces symmetry in such a way that the tightest lower bound is obtained.
#'
#' @references
#'
#' Keogh E and Ratanamahatana CA (2005). ``Exact indexing of dynamic time warping.'' *Knowledge and
#' information systems*, **7**(3), pp. 358-386.
#'
#' @examples
#'
#' # Sample data
#' data(uciCT)
#'
#' # Lower bound distance between two series
#' d.lbk <- lb_keogh(CharTraj[[1]], CharTraj[[2]], window.size = 20)$d
#'
#' # Corresponding true DTW distance
#' d.dtw <- dtw(CharTraj[[1]], CharTraj[[2]],
#'              window.type = "slantedband", window.size = 20)$distance
#'
#' d.lbk <= d.dtw
#'
#' # Calculating the LB between several time series using the 'proxy' package
#' # (notice how both argments must be lists)
#' D.lbk <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "LB_Keogh",
#'                      window.size = 20, norm = "L2")
#'
#' # Corresponding true DTW distance
#' # (see dtwclust-package description for an explanation of DTW2)
#' D.dtw <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "DTW2",
#'                      window.type = "slantedband", window.size = 20)
#'
#' D.lbk <= D.dtw
#'
lb_keogh <- function(x, y, window.size = NULL, norm = "L1",
                     lower.env = NULL, upper.env = NULL,
                     force.symmetry = FALSE, error.check = TRUE)
{
    norm <- match.arg(norm, c("L1", "L2"))

    if (error.check) check_consistency(x, "ts")

    if (is_multivariate(list(x)))
        stop("lb_keogh does not support multivariate series.")

    if (is.null(lower.env) || is.null(upper.env)) {
        if (error.check) check_consistency(y, "ts")

        if (is_multivariate(list(y)))
            stop("lb_keogh does not support multivariate series.")

        if (length(x) != length(y))
            stop("The series must have the same length")

        window.size <- check_consistency(window.size, "window")
    }

    ## NOTE: the 'window.size' definition varies betwen dtw/call_envelop and runmin/max
    if (is.null(lower.env) && is.null(upper.env)) {
        envelopes <- compute_envelop(y, window.size = window.size, error.check = FALSE)
        lower.env <- envelopes$lower
        upper.env <- envelopes$upper

    } else if (is.null(lower.env)) {
        lower.env <- caTools::runmin(y, window.size*2L + 1L)

    } else if (is.null(upper.env)) {
        upper.env <- caTools::runmax(y, window.size*2L + 1L)
    }

    if (length(lower.env) != length(x))
        stop("Length mismatch between 'x' and the lower envelop")

    if (length(upper.env) != length(x))
        stop("Length mismatch between 'x' and the upper envelop")

    D <- rep(0, length(x))

    ind1 <- x > upper.env
    D[ind1] <- x[ind1] - upper.env[ind1]
    ind2 <- x < lower.env
    D[ind2] <- lower.env[ind2] - x[ind2]

    d <- switch(EXPR = norm,
                L1 = sum(D),
                L2 = sqrt(sum(D^2)))

    if (force.symmetry) {
        d2 <- lb_keogh(x = y, y = x, window.size = window.size, norm = norm)

        if (d2$d > d) {
            d <- d2$d
            lower.env = d2$lower.env
            upper.env = d2$upper.env
        }
    }

    ## Finish
    list(d = d,
         upper.env = upper.env,
         lower.env = lower.env)
}

# ========================================================================================================
# Loop without using native 'proxy' looping (to avoid multiple calculations of the envelop)
# ========================================================================================================

lb_keogh_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1", ...,
                           force.symmetry = FALSE, pairwise = FALSE, error.check = TRUE)
{
    norm <- match.arg(norm, c("L1", "L2"))

    window.size <- check_consistency(window.size, "window")

    x <- any2list(x)

    if (error.check)
        check_consistency(x, "tslist")

    if (is.null(y)) {
        y <- x

    } else {
        y <- any2list(y)

        if (error.check) check_consistency(y, "tslist")
    }

    if (is_multivariate(x) || is_multivariate(y))
        stop("lb_keogh does not support multivariate series.")

    retclass <- "crossdist"

    envelops <- lapply(y, function(s) { compute_envelop(s, window.size, error.check = FALSE) })

    lower.env <- lapply(envelops, "[[", "lower")
    upper.env <- lapply(envelops, "[[", "upper")

    lower.env <- split_parallel(lower.env)
    upper.env <- split_parallel(upper.env)

    if (pairwise) {
        X <- split_parallel(x)

        validate_pairwise(X, lower.env)

        D <- foreach(x = X, lower.env = lower.env, upper.env = upper.env,
                     .packages = "dtwclust",
                     .combine = c,
                     .multicombine = TRUE) %op% {
                         mapply(upper.env, lower.env, x,
                                FUN = function(u, l, x) {
                                    lb_keogh(x,
                                             norm = norm,
                                             lower.env = l,
                                             upper.env = u,
                                             error.check = FALSE)$d
                                })
                     }

        retclass <- "pairdist"

    } else {
        D <- foreach(lower.env = lower.env, upper.env = upper.env,
                     .packages = "dtwclust",
                     .combine = cbind,
                     .multicombine = TRUE) %op% {
                         ret <- mapply(U = upper.env, L = lower.env,
                                       MoreArgs = list(x = x),
                                       SIMPLIFY = FALSE,
                                       FUN = function(U, L, x) {
                                           ## This will return one column of the distance matrix
                                           D <- sapply(x, u = U, l = L,
                                                       FUN = function(x, u, l) {
                                                           lb_keogh(x,
                                                                    norm = norm,
                                                                    lower.env = l,
                                                                    upper.env = u,
                                                                    error.check = FALSE)$d
                                                       })
                                           D
                                       })

                         do.call(cbind, ret)
                     }
    }

    if (force.symmetry && !pairwise) {
        if (nrow(D) != ncol(D)) {
            warning("Unable to force symmetry. Resulting distance matrix is not square.")

        } else {
            ind.tri <- lower.tri(D)

            new.low.tri.vals <- t(D)[ind.tri]
            indCorrect <- D[ind.tri] > new.low.tri.vals
            new.low.tri.vals[indCorrect] <- D[ind.tri][indCorrect]

            D[ind.tri] <- new.low.tri.vals
            D <- t(D)
            D[ind.tri] <- new.low.tri.vals
        }
    }

    class(D) <- retclass
    attr(D, "method") <- "LB_Keogh"

    D
}
