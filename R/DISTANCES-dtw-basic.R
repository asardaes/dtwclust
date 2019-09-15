#' Basic DTW distance
#'
#' This is a custom implementation of the DTW algorithm without all the functionality included in
#' [dtw::dtw()]. Because of that, it should be faster, while still supporting the most common
#' options.
#'
#' @export
#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#'
#' @param x,y Time series. Multivariate series must have time spanning the rows and variables
#'   spanning the columns.
#' @param window.size Size for slanted band window. `NULL` means no constraint.
#' @param norm Norm for the LCM calculation, "L1" for Manhattan or "L2" for (squared) Euclidean. See
#'   notes.
#' @param step.pattern Step pattern for DTW. Only `symmetric1` or `symmetric2` supported here. Note
#'   that these are *not* characters. See [dtw::stepPattern].
#' @param backtrack Also compute the warping path between series? See details.
#' @param normalize Should the distance be normalized? Only supported for `symmetric2`.
#' @param sqrt.dist Only relevant for `norm = "L2"`, see notes.
#' @param ... Currently ignored.
#' @template error-check
#'
#' @details
#'
#' If `backtrack` is `TRUE`, the mapping of indices between series is returned in a list.
#'
#' @template window
#'
#' @return The DTW distance. For `backtrack` `=` `TRUE`, a list with:
#'
#'   - `distance`: The DTW distance.
#'   - `index1`: `x` indices for the matched elements in the warping path.
#'   - `index2`: `y` indices for the matched elements in the warping path.
#'
#' @template proxy
#' @template symmetric
#' @section Proxy version:
#'
#'   In order for symmetry to apply here, the following must be true: no window constraint is used
#'   (`window.size` is `NULL`) or, if one is used, all series have the same length.
#'
#' @note
#'
#' The elements of the local cost matrix are calculated by using either Manhattan or squared
#' Euclidean distance. This is determined by the `norm` parameter. When the squared Euclidean
#' version is used, the square root of the resulting DTW distance is calculated at the end (as
#' defined in Ratanamahatana and Keogh 2004; Lemire 2009; see vignette references). This can be
#' avoided by passing `FALSE` in `sqrt.dist`.
#'
#' The DTW algorithm (and the functions that depend on it) might return different values in 32 bit
#' installations compared to 64 bit ones.
#'
#' An infinite distance value indicates that the constraints could not be fulfilled, probably due to
#' a too small `window.size` or a very large length difference between the series.
#'
#' @example man-examples/multivariate-dtw.R
#'
dtw_basic <- function(x, y, window.size = NULL, norm = "L1",
                      step.pattern = dtw::symmetric2, backtrack = FALSE,
                      normalize = FALSE, sqrt.dist = TRUE, ..., error.check = TRUE)
{
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }

    if (is.null(window.size))
        window.size <- -1L
    else
        window.size <- check_consistency(window.size, "window")

    if (NCOL(x) != NCOL(y)) stop("Multivariate series must have the same number of variables.")

    if (identical(step.pattern, dtw::symmetric1))
        step.pattern <- 1
    else if (identical(step.pattern, dtw::symmetric2))
        step.pattern <- 2
    else
        stop("step.pattern must be either symmetric1 or symmetric2 (without quotes)")

    norm <- match.arg(norm, c("L1", "L2"))
    norm <- switch(norm, "L1" = 1, "L2" = 2)
    backtrack <- isTRUE(backtrack)
    normalize <- isTRUE(normalize)
    sqrt.dist <- isTRUE(sqrt.dist)
    if (normalize && step.pattern == 1) stop("Unable to normalize with chosen step pattern.")

    if (backtrack)
        gcm <- matrix(0, NROW(x) + 1L, NROW(y) + 1L)
    else
        gcm <- matrix(0, 2L, NROW(y) + 1L)

    d <- .Call(C_dtw_basic, x, y, window.size,
               NROW(x), NROW(y), NCOL(x),
               norm, step.pattern, backtrack, normalize, sqrt.dist,
               gcm, PACKAGE = "dtwclust")

    if (backtrack) {
        d$index1 <- d$index1[d$path:1L]
        d$index2 <- d$index2[d$path:1L]
        d$path <- NULL
    }
    # return
    d
}

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#'
dtw_basic_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1",
                            step.pattern = dtw::symmetric2,
                            normalize = FALSE, sqrt.dist = TRUE, ...,
                            error.check = TRUE, pairwise = FALSE)
{
    x <- tslist(x)
    if (error.check) check_consistency(x, "vltslist")
    if (is.null(y)) {
        y <- x
        symmetric <- is.null(window.size) || !different_lengths(x)
    }
    else {
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
        symmetric <- FALSE
    }

    fill_type <- mat_type <- dim_names <- NULL # avoid warning about undefined globals
    eval(prepare_expr) # UTILS-expressions.R

    # adjust parameters for this distance
    if (is.null(window.size))
        window.size <- -1L
    else
        window.size <- check_consistency(window.size, "window")

    if (identical(step.pattern, dtw::symmetric1))
        step.pattern <- 1
    else if (identical(step.pattern, dtw::symmetric2))
        step.pattern <- 2
    else
        stop("step.pattern must be either symmetric1 or symmetric2 (without quotes)")

    normalize <- isTRUE(normalize)
    sqrt.dist <- isTRUE(sqrt.dist)

    if (normalize && step.pattern == 1) stop("Unable to normalize with chosen step pattern.")

    norm <- match.arg(norm, c("L1", "L2"))
    norm <- switch(norm, "L1" = 1, "L2" = 2)

    mv <- is_multivariate(c(x, y))
    backtrack <- FALSE

    # calculate distance matrix
    distance <- "DTW_BASIC" # read in C++, can't be temporary!
    distargs <- list(
        window.size = window.size,
        norm = norm,
        step.pattern = step.pattern,
        backtrack = backtrack,
        normalize = normalize,
        sqrt.dist = sqrt.dist
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
    attr(D, "method") <- "DTW_BASIC"
    # return
    D
}
