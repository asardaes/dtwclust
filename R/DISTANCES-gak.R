#' @importFrom stats median
#'
estimate_sigma <- function(x, y, within_proxy) {
    if (within_proxy)
        L <- c(sapply(x, NROW), sapply(y, NROW))
    else
        L <- c(NROW(x), NROW(y))

    pool <- unlist(c(x,y))
    rep <- stats::median(L)
    n <- round(0.5 * min(L))
    med1 <- sqrt(stats::median(L))
    med2 <- stats::median(replicate(rep, {
        xx <- sample(pool, n)
        yy <- sample(pool, n)
        l2norm(xx - yy) # UTILS-utils.R
    }))
    # return
    med1 * med2
}

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
#' @param error.check `r roxygen_error_check_param()`
#'
#' @details
#'
#' This function uses the Triangular Global Alignment Kernel (TGAK) described in Cuturi (2011). It
#' supports series of different length and multivariate series, so long as the ratio of the series'
#' lengths doesn't differ by more than 2 (or less than 0.5).
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
#' @section `r roxygen_proxy_section()`
#'
#' `r roxygen_proxy_symmetric()`
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
GAK <- function(x, y, ..., sigma = NULL, window.size = NULL, normalize = TRUE, error.check = TRUE)
{
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }

    # checks dimension consistency
    if (is_multivariate(list(x,y))) {
        x <- base::as.matrix(x)
        y <- base::as.matrix(y)
    }
    else {
        x <- as.numeric(x)
        y <- as.numeric(y)
    }

    if (is.null(window.size))
        window.size <- 0L
    else
        window.size <- check_consistency(window.size, "window")

    if (is.null(sigma))
        sigma <- estimate_sigma(x, y, FALSE)
    else if (sigma <= 0)
        stop("Parameter 'sigma' must be positive.")

    logs <- matrix(0, max(NROW(x), NROW(y)) + 1L, 3L)
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
    # return
    logGAK
}

#' @rdname GAK
#' @export
#'
gak <- GAK

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

gak_proxy <- function(x, y = NULL, ..., sigma = NULL, window.size = NULL, normalize = TRUE,
                      error.check = TRUE, pairwise = FALSE, .internal_ = FALSE)
{
    # normalization will be done manually to avoid multiple calculations of gak_x and gak_y
    if (!.internal_ && !normalize) { # nocov start
        warning("The proxy::dist version of GAK is always normalized.")
        normalize <- TRUE
    } # nocov end

    x <- tslist(x)
    if (error.check) check_consistency(x, "vltslist")

    if (is.null(y)) {
        symmetric <- normalize
        y <- x
    }
    else {
        symmetric <- FALSE
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
    }

    if (is.null(sigma))
        sigma <- estimate_sigma(x, y, TRUE)
    else if (sigma <= 0)
        stop("Parameter 'sigma' must be positive.")

    fill_type <- mat_type <- dim_names <- NULL # avoid warning about undefined globals
    eval(prepare_expr) # UTILS-expressions.R

    # adjust parameters for this distance
    if (normalize) {
        # calculation of normalization factors
        # x
        gak_x <- gak_proxy(x, x,
                           sigma = sigma, window.size = window.size,
                           error.check = FALSE, pairwise = TRUE, normalize = FALSE,
                           .internal_ = TRUE)
        # y
        if (symmetric)
            gak_y <- gak_x
        else
            gak_y <- gak_proxy(y, y,
                               sigma = sigma, window.size = window.size,
                               error.check = FALSE, pairwise = TRUE, normalize = FALSE,
                               .internal_ = TRUE)
    }

    mv <- is_multivariate(c(x,y))

    if (is.null(window.size))
        window.size <- 0L
    else
        window.size <- check_consistency(window.size, "window")

    # calculate distance matrix
    distance <- "GAK" # read in C++, can't be temporary!
    distargs <- list(
        sigma = sigma,
        window.size = window.size
    )
    num_threads <- get_nthreads()
    .Call(C_distmat_loop,
          D, x, y, distance, distargs, fill_type, mat_type, num_threads,
          PACKAGE = "dtwclust")

    # adjust D's attributes
    if (pairwise) {
        dim(D) <- NULL
        if (normalize) D <- 1 - exp(D - 0.5 * (gak_x + gak_y))
        class(D) <- "pairdist"
    }
    else if (symmetric) {
        dim(D) <- NULL

        # normalize
        j_upper <- length(x) - 1L
        i_lower <- 1L
        k <- 1L
        for (j in 1L:j_upper) {
            for (i in (j + i_lower):length(x)) {
                if (i != j) {
                    D[k] <- 1 - exp(D[k] - (gak_x[i] + gak_x[j]) / 2)
                }
                k <- k + 1L
            }
        }

        class(D) <- "dist"
        attr(D, "Size") <- length(x)
        attr(D, "Labels") <- names(x)
    }
    else {
        if (normalize) D <- 1 - exp(D - outer(gak_x, gak_y, function(x, y) { (x + y) / 2 }))
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }


    attr(D, "method") <- "GAK"
    attr(D, "sigma") <- sigma
    # return
    D
}

# ==================================================================================================
# Wrapper for proxy::simil
# ==================================================================================================

gak_simil <- function(x, y = NULL, ..., normalize = FALSE) {
    if (normalize) warning("The proxy::simil version of GAK cannot be normalized.") # nocov
    gak_proxy(x = x, y = y, ..., normalize = FALSE, .internal_ = TRUE)
}
