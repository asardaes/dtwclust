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
#' @template error-check
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
#' @return
#'
#' The logarithm of the GAK if `normalize = FALSE`, otherwise 1 minus the normalized GAK. The value
#' of `sigma` is assigned as an attribute of the result.
#'
#' @template proxy
#' @template symmetric
#'
#' @note
#'
#' The estimation of `sigma` does *not* depend on `window.size`.
#'
#' If `normalize` is set to `FALSE`, the returned value is **not** a distance, rather a similarity.
#' The [proxy::dist()] version is thus always normalized. Use [proxy::simil()] with `method` set to
#' "uGAK" if you want the unnormalized similarities.
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
#' # Unnormalized similarities
#' proxy::simil(CharTraj[1L:5L], method = "ugak")
#'
GAK <- function(x, y, ..., sigma = NULL, window.size = NULL, normalize = TRUE,
                logs = NULL, error.check = TRUE)
{
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

    if (is.null(window.size))
        window.size <- 0L
    else
        window.size <- check_consistency(window.size, "window")

    if (is.null(sigma)) {
        med1 <- sqrt(median(c(NROW(x), NROW(y))))
        n <- round(0.5 * min(NROW(x), NROW(y)))

        if (is.null(dim(x)))
            xx <- sample(x, n)
        else
            xx <- sapply(1L:ncol(x), function(idc) { sample(x[, idc], n) })

        if (is.null(dim(y)))
            yy <- sample(y, n)
        else
            yy <- sapply(1L:ncol(y), function(idc) { sample(y[, idc], n) })

        med2 <- median(replicate(100L, lnorm(xx - yy, 2)))
        sigma <- med1 * med2

    } else if (sigma <= 0) stop("Parameter 'sigma' must be positive.")

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

    ## return
    logGAK
}

#' @rdname GAK
#' @export
#'
gak <- GAK

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

GAK_proxy <- function(x, y = NULL, ..., sigma = NULL, window.size = NULL, normalize = TRUE,
                      logs = NULL, error.check = TRUE, pairwise = FALSE, .internal_ = FALSE)
{
    ## normalization will be done manually to avoid multiple calculations of gak_x and gak_y
    if (!.internal_ && !normalize) { # nocov start
        warning("The proxy::dist version of GAK is always normalized.")
        normalize <- TRUE
    } # nocov end

    x <- tslist(x)
    if (error.check) check_consistency(x, "vltslist")

    if (is.null(y)) {
        symmetric <- normalize
        y <- x

    } else {
        symmetric <- FALSE
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
    }

    if (is.null(sigma)) {
        L <- c(sapply(x, NROW), sapply(y, NROW))
        n <- round(0.5 * min(L))
        med1 <- sqrt(median(L))

        med2 <- median(replicate(length(x) + length(y), {
            x <- sample(x, 1L)[[1L]]
            y <- sample(y, 1L)[[1L]]

            if (is.null(dim(x)))
                xx <- sample(x, n)
            else
                xx <- sapply(1L:ncol(x), function(idc) { sample(x[, idc], n) })

            if (is.null(dim(y)))
                yy <- sample(y, n)
            else
                yy <- sapply(1L:ncol(y), function(idc) { sample(y[, idc], n) })

            lnorm(xx - yy, 2)
        }))

        sigma <- med1 * med2

    } else if (sigma <= 0) stop("Parameter 'sigma' must be positive.")

    ## parallel chunks are made column-wise, so flip x and y if necessary
    flip <- NULL
    num_workers <- foreach::getDoParWorkers()

    if (!pairwise && !symmetric && length(y) < num_workers && length(x) >= num_workers) {
        flip <- y
        y <- x
        x <- flip
    }

    ## pre-allocate logs
    if (is.null(logs)) logs <- matrix(0, max(sapply(x, NROW), sapply(y, NROW)) + 1L, 3L)
    pairwise <- isTRUE(pairwise)
    dim_out <- c(length(x), length(y))
    dim_names <- list(names(x), names(y))
    D <- allocate_distmat(length(x), length(y), pairwise, symmetric) ## utils.R

    if (normalize) {
        ## calculation of normalization factors
        # x
        gak_x <- GAK_proxy(x, x,
                           sigma = sigma, window.size = window.size, logs = logs,
                           error.check = FALSE, pairwise = TRUE, normalize = FALSE,
                           .internal_ = TRUE)
        # y
        if (symmetric)
            gak_y <- gak_x
        else
            gak_y <- GAK_proxy(y, y,
                               sigma = sigma, window.size = window.size, logs = logs,
                               error.check = FALSE, pairwise = TRUE, normalize = FALSE,
                               .internal_ = TRUE)
    }

    ## Wrap as needed for foreach
    if (pairwise) {
        x <- split_parallel(x)
        y <- split_parallel(y)
        validate_pairwise(x, y)
        endpoints <- attr(x, "endpoints")

    } else if (symmetric) {
        endpoints <- symmetric_loop_endpoints(length(x)) ## utils.R
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- x

    } else {
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- split_parallel(y)
        endpoints <- attr(y, "endpoints")
    }

    if (bigmemory::is.big.matrix(D)) {
        D_desc <- bigmemory::describe(D)
        noexport <- c("D", "gak_x", "gak_y")
        packages <- c("dtwclust", "bigmemory")

    } else {
        D_desc <- NULL
        noexport <- c("gak_x", "gak_y")
        packages <- c("dtwclust")
    }

    ## Calculate distance matrix
    foreach(x = x, y = y, endpoints = endpoints,
            .combine = c,
            .multicombine = TRUE,
            .packages = packages,
            .export = c("gak_loop"),
            .noexport = noexport) %op% {
                bigmat <- !is.null(D_desc)
                d <- if (bigmat) bigmemory::attach.big.matrix(D_desc)@address else D
                do.call(gak_loop,
                        list(d = d,
                             x = x,
                             y = y,
                             symmetric = symmetric,
                             pairwise = pairwise,
                             endpoints = endpoints,
                             bigmat = bigmat,
                             window.size = window.size,
                             sigma = sigma,
                             logs = logs),
                        TRUE)
            }

    D <- D[,]
    if (pairwise) {
        if (normalize) D <- 1 - exp(D - 0.5 * (gak_x + gak_y))
        class(D) <- "pairdist"

    } else {
        if (is.null(dim(D))) dim(D) <- dim_out
        if (normalize) D <- 1 - exp(D - outer(gak_x, gak_y, function(x, y) { (x + y) / 2 }))
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }

    if (!pairwise && symmetric && normalize) diag(D) <- 0
    if (!is.null(flip)) D <- t(D)
    attr(D, "method") <- "GAK"
    attr(D, "sigma") <- sigma
    ## return
    D
}

# ==================================================================================================
# Wrapper for C++
# ==================================================================================================

gak_loop <- function(d, x, y, symmetric, pairwise, endpoints, bigmat, ...,
                     window.size, sigma, logs)
{
    ## check dimension consistency
    mv <- is_multivariate(c(x,y))

    nr <- max(sapply(x, NROW), sapply(y, NROW)) + 1L
    if (is.null(logs))
        logs <- matrix(0, nr, 3L)
    else if (!is.matrix(logs) || nrow(logs) < nr || ncol(logs) < 3L)
        stop("GAK: Dimension inconsistency in 'logs'")
    else if (storage.mode(logs) != "double")
        stop("GAK: If provided, 'logs' must have 'double' storage mode.")

    if (is.null(window.size))
        window.size <- 0L
    else
        window.size <- check_consistency(window.size, "window")

    distargs <- list(window.size = window.size,
                     sigma = sigma,
                     logs = logs)

    .Call(C_gak_loop,
          d, x, y, symmetric, pairwise, bigmat, mv, distargs, endpoints,
          PACKAGE = "dtwclust")
}

# ==================================================================================================
# Wrapper for proxy::simil
# ==================================================================================================

GAK_simil <- function(x, y = NULL, ..., normalize = FALSE) {
    if (normalize) warning("The proxy::simil version of GAK cannot be normalized.")
    GAK_proxy(x = x, y = y, ..., normalize = FALSE, .internal_ = TRUE)
}
