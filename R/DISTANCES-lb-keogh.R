#' Keogh's DTW lower bound
#'
#' This function calculates a lower bound (LB) on the Dynamic Time Warp (DTW) distance between two
#' time series. It uses a Sakoe-Chiba constraint.
#'
#' @export
#' @include DISTANCES-lb-improved.R
#'
#' @inheritParams lb_improved
#' @inherit lb_improved details
#' @inheritSection lb_improved Note
#'
#' @section `r roxygen_proxy_section()`
#'
#' @return A list with:
#'
#' - `d`: The lower bound of the DTW distance.
#' - `upper.env`: The time series of `y`'s upper envelope.
#' - `lower.env`: The time series of `y`'s lower envelope.
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
#'              window.type = "sakoechiba", window.size = 20)$distance
#'
#' d.lbk <= d.dtw
#'
#' # Calculating the LB between several time series using the 'proxy' package
#' # (notice how both argments must be lists)
#' D.lbk <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "LB_Keogh",
#'                      window.size = 20, norm = "L2")
#'
#' # Corresponding true DTW distance
#' D.dtw <- proxy::dist(CharTraj[1], CharTraj[2:5], method = "dtw_basic",
#'                      norm = "L2", window.size = 20)
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
    }
    else {
        check_consistency(lower.env, "ts")
        check_consistency(upper.env, "ts")
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
            lower.env <- d2$lower.env
            upper.env <- d2$upper.env
        }
    }
    # return
    list(d = d, upper.env = upper.env, lower.env = lower.env)
}

# ==================================================================================================
# Loop without using native 'proxy' looping (to avoid multiple calculations of the envelope)
# ==================================================================================================

lb_keogh_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1", ...,
                           force.symmetry = FALSE, pairwise = FALSE, error.check = TRUE)
{
    x <- tslist(x)
    if (is.null(y))
        y <- x
    else
        y <- tslist(y)
    if (length(x) == 0L || length(y) == 0L) stop("Empty list received in x or y.") # nocov start
    if (error.check) check_consistency(c(x,y), "tslist")
    if (is_multivariate(c(x,y))) stop("lb_keogh does not support multivariate series.") # nocov end

    fill_type <- mat_type <- dim_names <- NULL # avoid warning about undefined globals
    symmetric <- FALSE
    eval(prepare_expr) # UTILS-expressions.R

    # adjust parameters for this distance
    norm <- match.arg(norm, c("L1", "L2"))
    window.size <- check_consistency(window.size, "window")
    envelopes <- lapply(y, function(s) { compute_envelope(s, window.size, error.check = FALSE) })
    lower.env <- lapply(envelopes, "[[", "lower")
    upper.env <- lapply(envelopes, "[[", "upper")

    # calculate distance matrix
    distance <- "LBK" # read in C++, can't be temporary!
    distargs <- list(
        p = switch(norm, "L1" = 1L, "L2" = 2L),
        len = length(x[[1L]]),
        lower.env = lower.env,
        upper.env = upper.env
    )
    num_threads <- get_nthreads()
    .Call(C_distmat_loop,
          D, x, y, distance, distargs, fill_type, mat_type, num_threads,
          PACKAGE = "dtwclust")

    # adjust D's attributes
    if (pairwise) {
        dim(D) <- NULL
        class(D) <- "pairdist"
    }
    else {
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }
    if (force.symmetry && !pairwise) {
        if (nrow(D) != ncol(D))
            warning("Unable to force symmetry. Resulting distance matrix is not square.") # nocov
        else
            .Call(C_force_lb_symmetry, D, PACKAGE = "dtwclust")
    }

    attr(D, "method") <- "LB_Keogh"
    D
}
