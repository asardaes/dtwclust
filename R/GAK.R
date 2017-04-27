#' Fast global alignment kernels
#'
#' Distance based on (triangular) global alignment kernels.
#'
#' @export
#'
#' @param x,y Time series. A multivariate series should have time spanning the rows and variables
#'   spanning the columns.
#' @param ... Currently ignored.
#' @param sigma Parameter for the Gaussian kernel's width. See details for the interpretation of
#'   `NULL`.
#' @param window.size Parameterization of the constraining band (*T* in Cuturi (2011)). See details.
#' @param normalize Normalize the result by considering diagonal terms.
#' @param logs Optionally, a matrix with `max(NROW(x), NROW(y)) + 1` rows and 3 columns to use for
#'   the logarithm calculations. Used internally for memory optimization. If provided, it **will**
#'   be modified *in place* by `C` code, except in the parallel version in [proxy::dist()] which
#'   ignores it for thread-safe reasons.
#' @param error.check Check data inconsistencies?
#'
#' @details
#'
#' This function uses the Triangular Global Alignment Kernel (TGAK) described in Cuturi (2011). It
#' supports series of different length and multivariate series, so long as the ratio of the series'
#' lengths don't differ by more than 2 (or less than 0.5).
#'
#' The `window.size` parameter is similar to the one used in DTW, so `NULL` signifies no constraint,
#' and its value should be greater than 1 if used with series of different length.
#'
#' The Gaussian kernel is parameterized by `sigma`. Providing `NULL` means that the value will be
#' estimated by using the strategy mentioned in Cuturi (2011) with a constant of 1. This estimation
#' is subject to **randomness**, so consider estimating the value once and re-using it (the estimate
#' is returned as an attribute of the result). See the examples.
#'
#' For more information, refer to the package vignette and the referenced article.
#'
#' @return The logarithm of the GAK if `normalize = FALSE`, otherwise 1 minus the normalized GAK.
#'   The value of `sigma` is assigned as an attribute of the result.
#'
#' @note
#'
#' If `normalize` is set to `FALSE`, the returned value is **not** a distance, rather a similarity.
#' The [proxy::dist()] version is thus always normalized.
#'
#' A constrained unnormalized calculation (i.e. with `window.size > 0` and `normalize = FALSE`) will
#' return negative infinity if `abs(NROW(x)` `-` `NROW(y))` `>` `window.size`. Since the function
#' won't perform calculations in that case, it might be faster, but if this behavior is not desired,
#' consider reinterpolating the time series (see [reinterpolate()]) or increasing the window size.
#'
#' @references
#'
#' Cuturi, M. (2011). Fast global alignment kernels. In *Proceedings of the 28th international
#' conference on machine learning (ICML-11)* (pp. 929-936).
#'
#' @examples
#'
#' \dontrun{
#' data(uciCT)
#'
#' set.seed(832)
#' GAKd <- proxy::dist(zscore(CharTraj), method = "gak",
#'                     pairwise = TRUE, window.size = 18L)
#'
#' # Obtained estimate of sigma
#' sigma <- attr(GAKd, "sigma")
#'
#' # Use value for clustering
#' tsclust(CharTraj, k = 20L,
#'         distance = "gak", centroid = "shape",
#'         trace = TRUE,
#'         args = tsclust_args(dist = list(sigma = sigma,
#'                                         window.size = 18L)))
#' }
#'
GAK <- function(x, y, ..., sigma = NULL, window.size = NULL, normalize = TRUE,
                logs = NULL, error.check = TRUE)
{
    x <- cbind(x)
    y <- cbind(y)

    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }

    ## check dimension consistency
    is_multivariate(list(x,y))

    if (is.null(logs))
        logs <- matrix(0, max(NROW(x), NROW(y)) + 1L, 3L)
    else if (!is.matrix(logs) || nrow(logs) < (max(NROW(x), NROW(y)) + 1L) || ncol(logs) < 3L)
        stop("GAK: Dimension inconsistency in 'logs'")
    else if (storage.mode(logs) != "double")
        stop("GAK: If provided, 'logs' must have 'double' storage mode.")

    if (is.null(sigma)) {
        n <- round(0.5 * min(NROW(x), NROW(y)))
        med1 <- sqrt(median(c(NROW(x), NROW(y))))

        xx <- sapply(1L:NCOL(x), function(idc) sample(x[ , idc], n))
        yy <- sapply(1L:NCOL(y), function(idc) sample(y[ , idc], n))

        med2 <- median(replicate(100L, lnorm(xx - yy, 2)))

        sigma <- med1 * med2

    } else if (sigma <= 0)
        stop("Parameter 'sigma' must be positive.")

    if (is.null(window.size))
        window.size <- 0L
    else
        window.size <- check_consistency(window.size, "window")

    logGAK <- .Call(C_logGAK, x, y,
                    NROW(x), NROW(y), NCOL(x),
                    sigma, window.size, logs,
                    PACKAGE = "dtwclust")

    if (normalize) {
        gak_x <- .Call(C_logGAK, x, x,
                       NROW(x), NROW(x), NCOL(x),
                       sigma, window.size, logs,
                       PACKAGE = "dtwclust")

        gak_y <- .Call(C_logGAK, y, y,
                       NROW(y), NROW(y), NCOL(y),
                       sigma, window.size, logs,
                       PACKAGE = "dtwclust")

        logGAK <- 1 - exp(logGAK - 0.5 * (gak_x + gak_y))
    }

    attr(logGAK, "sigma") <- sigma

    logGAK
}

GAK_proxy <- function(x, y = NULL, ..., sigma = NULL, normalize = TRUE, logs = NULL,
                      error.check = TRUE, pairwise = FALSE)
{
    if (!normalize)
        warning("The proxy version of GAK is always normalized.")

    x <- any2list(x)
    if (error.check) check_consistency(x, "vltslist")

    dots <- list(...)
    dots$error.check <- FALSE

    ## normalization will be done manually to avoid multiple calculations of gak_x and gak_y
    dots$normalize <- FALSE

    if (is.null(y)) {
        symmetric <- TRUE
        y <- x

    } else {
        symmetric <- FALSE
        y <- any2list(y)
        if (error.check) check_consistency(y, "vltslist")
    }

    if (is.null(sigma)) {
        L <- c(sapply(x, NROW), sapply(y, NROW))
        n <- round(0.5 * min(L))

        med1 <- sqrt(median(L))

        med2 <- median(replicate(length(x) + length(y), {
            x <- cbind(sample(x, 1L)[[1L]])
            y <- cbind(sample(y, 1L)[[1L]])

            xx <- sapply(1L:NCOL(x), function(idc) sample(x[ , idc], n))
            yy <- sapply(1L:NCOL(y), function(idc) sample(y[ , idc], n))

            lnorm(xx - yy, 2)
        }))

        sigma <- med1 * med2

    } else if (sigma <= 0)
        stop("Parameter 'sigma' must be positive.")

    dots$sigma <- sigma

    ## parallel chunks are made column-wise, so flip x and y if necessary
    flip <- NULL
    num_workers <- foreach::getDoParWorkers()

    if (!pairwise && !symmetric && length(y) < num_workers && length(x) >= num_workers) {
        flip <- y
        y <- x
        x <- flip
    }

    retclass <- "crossdist"

    ## to pre-allocate LOGS
    L <- max(sapply(x, NROW), sapply(y, NROW)) + 1L

    X <- split_parallel(x)
    Y <- split_parallel(y)

    ## calculation of normalization factors
    # x
    LOGS <- allocate_matrices(logs, nrow = L, ncol = 3L, target.size = length(X))

    gak_x <- foreach(x = X, logs = LOGS,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         sapply(x, function(xx) {
                             do.call("GAK",
                                     enlist(x = xx,
                                            y = xx,
                                            logs = logs,
                                            dots = dots))
                         })
                     }

    # y
    LOGS <- allocate_matrices(logs, nrow = L, ncol = 3L, target.size = length(Y))

    if (symmetric) {
        gak_y <- gak_x

    } else {
        gak_y <- foreach(y = Y, logs = LOGS,
                         .combine = c,
                         .multicombine = TRUE,
                         .packages = "dtwclust",
                         .export = "enlist") %op% {
                             sapply(y, function(yy) {
                                 do.call("GAK",
                                         enlist(x = yy,
                                                y = yy,
                                                logs = logs,
                                                dots = dots))
                             })
                         }
    }

    ## Calculate distance matrix
    if (pairwise) {
        validate_pairwise(X, Y)

        D <- foreach(x = X, y = Y, logs = LOGS,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         mapply(x, y, FUN = function(x, y) {
                             do.call("GAK",
                                     enlist(x = x,
                                            y = y,
                                            logs = logs,
                                            dots = dots))
                         })
                     }

        ## normalize
        D <- 1 - exp(D - 0.5 * (gak_x + gak_y))

        names(D) <- NULL
        retclass <- "pairdist"

    } else if (symmetric) {
        pairs <- call_pairs(length(x), lower = FALSE)
        pairs <- split_parallel(pairs, 1L)

        LOGS <- allocate_matrices(logs, nrow = L, ncol = 3L, target.size = length(pairs))

        d <- foreach(pairs = pairs, logs = LOGS,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         mapply(x[pairs[ , 1L]], x[pairs[ , 2L]],
                                SIMPLIFY = TRUE,
                                FUN = function(xx, yy) {
                                    do.call("GAK",
                                            enlist(x = xx,
                                                   y = yy,
                                                   logs = logs,
                                                   dots = dots))
                                })
                     }

        rm("pairs")

        D <- matrix(0, nrow = length(x), ncol = length(x))
        D[upper.tri(D)] <- d
        D <- t(D)
        D[upper.tri(D)] <- d

        ## normalize
        D <- 1 - exp(D - outer(gak_x, gak_y, function(x, y) { (x + y) / 2 }))
        diag(D) <- 0

        attr(D, "dimnames") <- list(names(x), names(x))

    } else {
        D <- foreach(y = Y, logs = LOGS,
                     .combine = cbind,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         ret <- lapply(y, x = x, FUN = function(y, x) {
                             sapply(x, y = y, FUN = function(x, y) {
                                 do.call("GAK",
                                         enlist(x = x,
                                                y = y,
                                                logs = logs,
                                                dots = dots))
                             })
                         })

                         do.call(cbind, ret)
                     }

        ## normalize
        D <- 1 - exp(D - outer(gak_x, gak_y, function(x, y) { (x + y) / 2 }))
    }

    if (!is.null(flip)) D <- t(D)
    class(D) <- retclass
    attr(D, "method") <- "GAK"
    attr(D, "sigma") <- sigma

    D
}
