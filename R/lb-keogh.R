#' Keogh's DTW lower bound
#'
#' This function calculates a lower bound (LB) on the Dynamic Time Warp (DTW) distance between two
#' time series. It uses a Sakoe-Chiba constraint.
#'
#' @export
#' @include lb-improved.R
#'
#' @inheritParams lb_improved
#' @inherit lb_improved details
#' @inheritSection lb_improved Note
#'
#' @return A list with:
#'
#'   - `d`: The lower bound of the DTW distance.
#'   - `upper.env`: The time series of `y`'s upper envelope.
#'   - `lower.env`: The time series of `y`'s lower envelope.
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
    if (is_multivariate(list(x))) stop("lb_keogh does not support multivariate series.")

    if (is.null(lower.env) || is.null(upper.env)) {
        if (is_multivariate(list(y))) stop("lb_keogh does not support multivariate series.")
        if (length(x) != length(y)) stop("The series must have the same length")
        if (error.check) check_consistency(y, "ts")
        window.size <- check_consistency(window.size, "window")

        envelopes <- compute_envelope(y, window.size = window.size, error.check = FALSE)
        lower.env <- envelopes$lower
        upper.env <- envelopes$upper

    } else {
        if (length(lower.env) != length(x))
            stop("Length mismatch between 'x' and the lower envelope")
        if (length(upper.env) != length(x))
            stop("Length mismatch between 'x' and the upper envelope")
    }

    p <- switch(norm, L1 = 1L, L2 = 2L)
    d <- .Call(C_lbk, x, p, lower.env, upper.env, PACKAGE = "dtwclust")

    if (force.symmetry) {
        d2 <- lb_keogh(x = y, y = x, window.size = window.size, norm = norm, error.check = FALSE)
        if (d2$d > d) {
            d <- d2$d
            lower.env = d2$lower.env
            upper.env = d2$upper.env
        }
    }

    ## return
    list(d = d, upper.env = upper.env, lower.env = lower.env)
}

# ==================================================================================================
# Loop without using native 'proxy' looping (to avoid multiple calculations of the envelope)
# ==================================================================================================

lb_keogh_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1", ...,
                           force.symmetry = FALSE, pairwise = FALSE, error.check = TRUE)
{
    norm <- match.arg(norm, c("L1", "L2"))
    window.size <- check_consistency(window.size, "window")
    x <- any2list(x)
    if (error.check) check_consistency(x, "tslist")

    if (is.null(y)) {
        y <- x

    } else {
        y <- any2list(y)
        if (error.check) check_consistency(y, "tslist")
    }

    if (is_multivariate(x) || is_multivariate(y))
        stop("lb_keogh does not support multivariate series.")

    envelopes <- lapply(y, function(s) { compute_envelope(s, window.size, error.check = FALSE) })
    lower.env <- lapply(envelopes, "[[", "lower")
    upper.env <- lapply(envelopes, "[[", "upper")
    lower.env <- split_parallel(lower.env)
    upper.env <- split_parallel(upper.env)


    if (pairwise) {
        X <- split_parallel(x)
        validate_pairwise(X, lower.env)
        retclass <- "pairdist"

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

    } else {
        retclass <- "crossdist"

        D <- foreach(lower.env = lower.env, upper.env = upper.env,
                     .packages = "dtwclust",
                     .combine = cbind,
                     .multicombine = TRUE) %op% {
                         d <- mapply(U = upper.env, L = lower.env,
                                     MoreArgs = list(x = x),
                                     SIMPLIFY = FALSE,
                                     FUN = function(U, L, x) {
                                         ## This will return one column of the distance matrix
                                         sapply(x, u = U, l = L,
                                                FUN = function(x, u, l) {
                                                    lb_keogh(x,
                                                             norm = norm,
                                                             lower.env = l,
                                                             upper.env = u,
                                                             error.check = FALSE)$d
                                                })
                                     })

                         do.call(cbind, d)
                     }
    }

    if (force.symmetry && !pairwise) {
        if (nrow(D) != ncol(D)) {
            warning("Unable to force symmetry. Resulting distance matrix is not square.")

        } else {
            ind_tri <- lower.tri(D)
            new_low_tri_vals <- t(D)[ind_tri]
            ind_correct <- D[ind_tri] > new_low_tri_vals
            new_low_tri_vals[ind_correct] <- D[ind_tri][ind_correct]
            D[ind_tri] <- new_low_tri_vals
            D <- t(D)
            D[ind_tri] <- new_low_tri_vals
        }
    }

    class(D) <- retclass
    attr(D, "method") <- "LB_Keogh"

    ## return
    D
}
