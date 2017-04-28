#' Basic DTW distance
#'
#' This is a custom implementation of the DTW algorithm without all the functionality included in
#' [dtw::dtw()]. Because of that, it should be slightly faster, while still supporting the most
#' common options.
#'
#' @export
#'
#' @param x,y Time series. Multivariate series must have time spanning the rows and variables
#'   spanning the columns.
#' @param window.size Size for slanted band window. `NULL` means no constraint.
#' @param norm Norm for the DTW calculation, "L1" for Manhattan or "L2" for Euclidean.
#' @param step.pattern Step pattern for DTW. Only `symmetric1` or `symmetric2` supported here. Note
#'   that these are *not* characters. See [dtw::stepPattern].
#' @param backtrack Also compute the warping path between series? See details.
#' @param normalize Should the distance be normalized? Only supported for `symmetric2`.
#' @param ... Currently ignored.
#' @param gcm Optionally, a matrix to use for the global cost matrix calculations. It should have
#'   `NROW(y)+1` columns and `NROW(x)+1` rows for `backtrack = TRUE` **or** `2` rows for `backtrack
#'   = FALSE`. Used internally for memory optimization. If provided, it **will** be modified *in
#'   place* by `C` code, except in the parallel version in [proxy::dist()] which ignores it for
#'   thread-safe reasons.
#' @param error.check Check data inconsistencies?
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
#' @note
#'
#' The DTW algorithm (and the functions that depend on it) might return different values in 32 bit
#' installations compared to 64 bit ones.
#'
dtw_basic <- function(x, y, window.size = NULL, norm = "L1",
                      step.pattern = symmetric2, backtrack = FALSE,
                      normalize = FALSE, ..., gcm = NULL, error.check = TRUE)
{
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }

    backtrack <- isTRUE(backtrack)

    if (NCOL(x) != NCOL(y))
        stop("Multivariate series must have the same number of variables.")

    if (is.null(window.size))
        window.size <- -1L
    else
        window.size <- check_consistency(window.size, "window")

    norm <- match.arg(norm, c("L1", "L2"))
    norm <- switch(norm, "L1" = 1, "L2" = 2)

    if (identical(step.pattern, symmetric1))
        step.pattern <- 1
    else if (identical(step.pattern, symmetric2))
        step.pattern <- 2
    else
        stop("step.pattern must be either symmetric1 or symmetric2")

    if (normalize && step.pattern == 1)
        stop("Unable to normalize with chosen step pattern.")

    if (backtrack) {
        if (is.null(gcm))
            gcm <- matrix(0, NROW(x) + 1L, NROW(y) + 1L)
        else if (!is.matrix(gcm) || nrow(gcm) < (NROW(x) + 1L) || ncol(gcm) < (NROW(y) + 1L))
            stop("dtw_basic: Dimension inconsistency in 'gcm'")

    } else {
        if (is.null(gcm))
            gcm <- matrix(0, 2L, NROW(y) + 1L)
        else if (!is.matrix(gcm) || nrow(gcm) < 2L || ncol(gcm) < (NROW(y) + 1L))
            stop("dtw_basic: Dimension inconsistency in 'gcm'")
    }

    if (storage.mode(gcm) != "double")
        stop("dtw_basic: If provided, 'gcm' must have 'double' storage mode.")

    d <- .Call(C_dtw_basic, x, y, window.size,
               NROW(x), NROW(y), NCOL(x),
               norm, step.pattern, backtrack,
               gcm, PACKAGE = "dtwclust")

    if (normalize) {
        if (backtrack)
            d$distance <- d$distance / (NROW(x) + NROW(y))
        else
            d <- d / (NROW(x) + NROW(y))
    }

    if (backtrack) {
        d$index1 <- d$index1[d$path:1L]
        d$index2 <- d$index2[d$path:1L]
        d$path <- NULL
    }

    d
}

dtw_basic_proxy <- function(x, y = NULL, ..., gcm = NULL, error.check = TRUE, pairwise = FALSE) {
    x <- any2list(x)
    if (error.check) check_consistency(x, "vltslist")

    dots <- list(...)
    dots$backtrack <- FALSE
    dots$error.check <- FALSE

    if (is.null(y)) {
        y <- x
        symmetric <- is.null(dots$window.size) || !different_lengths(x)

    } else {
        y <- any2list(y)
        if (error.check) check_consistency(y, "vltslist")

        symmetric <- FALSE
    }

    retclass <- "crossdist"

    ## to pre-allocate GCMs
    nc <- max(sapply(y, NROW)) + 1L

    ## Calculate distance matrix
    if (pairwise) {
        X <- split_parallel(x)
        Y <- split_parallel(y)

        validate_pairwise(X, Y)

        GCM <- allocate_matrices(gcm, nrow = 2L, ncol = nc, target.size = length(X))

        D <- foreach(x = X, y = Y, gcm = GCM,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         mapply(x, y, FUN = function(x, y) {
                             do.call("dtw_basic",
                                     enlist(x = x,
                                            y = y,
                                            gcm = gcm,
                                            dots = dots))
                         })
                     }

        names(D) <- NULL
        retclass <- "pairdist"

    } else if (symmetric) {
        pairs <- call_pairs(length(x), lower = FALSE)
        pairs <- split_parallel(pairs, 1L)

        GCM <- allocate_matrices(gcm, nrow = 2L, ncol = nc, target.size = length(pairs))

        dots$pairwise <- TRUE

        d <- foreach(pairs = pairs, gcm = GCM,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         do.call(proxy::dist,
                                 enlist(x = x[pairs[ , 1L]],
                                        y = x[pairs[ , 2L]],
                                        method = "dtw_basic",
                                        gcm = gcm,
                                        dots = dots))
                     }

        rm("pairs")

        D <- matrix(0, nrow = length(x), ncol = length(x))
        D[upper.tri(D)] <- d
        D <- t(D)
        D[upper.tri(D)] <- d

        attr(D, "dimnames") <- list(names(x), names(x))

    } else {
        Y <- split_parallel(y)

        GCM <- allocate_matrices(gcm, nrow = 2L, ncol = nc, target.size = length(Y))

        D <- foreach(y = Y, gcm = GCM,
                     .combine = cbind,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         ret <- lapply(y, x = x, FUN = function(y, x) {
                             sapply(x, y = y, FUN = function(x, y) {
                                 do.call("dtw_basic",
                                         enlist(x = x,
                                                y = y,
                                                gcm = gcm,
                                                dots = dots))
                             })
                         })

                         do.call(cbind, ret)
                     }
    }

    class(D) <- retclass
    attr(D, "method") <- "DTW_BASIC"

    D
}
