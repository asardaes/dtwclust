#' Fast global alignment kernels
#'
#' Distance based on (triangular) global alignment kernels.
#'
#' @param x,y Time series. A multivariate series should have time spanning the rows and variables spanning
#' the columns.
#' @param ... Currently ignored.
#' @param sigma Parameter for the Gaussian kernel's width. See details for the interpretation of \code{NULL}.
#' @param window.size Parameterization of the constraining band (\emph{T} in Cuturi (2011)). See details.
#' @param normalize Normalize the result by considering diagonal terms.
#'
#' @export
#'
#' @details
#'
#' This function uses the Triangular Global Alignment Kernel (TGAK) described in Cuturi (2011). It supports
#' series of different length and multivariate series, so long as the ratio of the series' lengths don't differ
#' by more than 2 (or less than 0.5).
#'
#' The \code{window.size} parameter is similar to the one used in DTW, so \code{NULL} signifies no constraint,
#' and its value should be greater than 1 for series of different length.
#'
#' The Gaussian kernel is parameterized by \code{sigma}. Providing \code{NULL} means that the value will be
#' estimated by using the strategy mentioned in Cuturi (2011) with a constant of 1. This estimation is subject
#' to \strong{randomness}, so consider estimating the value once and re-using it (the estimate is returned as
#' an attribute of the result). See the examples.
#'
#' For more information, refer to the package vignette and the referenced article.
#'
#' @return The logarithm of the GAK if \code{normalize = FALSE}, otherwise 1 minus the normalized GAK.
#' The value of \code{sigma} is assigned as an attribute of the result.
#'
#' @note
#'
#' If \code{normalize} is set to \code{FALSE}, the returned value is \strong{not} a distance, rather a similarity.
#' The \code{proxy::}\code{\link[proxy]{dist}} version is thus always normalized.
#'
#' A constrained unnormalized calculation (i.e. with \code{window.size > 0} and \code{normalize = FALSE}) will
#' return negative infinity if \code{abs(NROW(x)} \code{-} \code{NROW(y))} \code{>} \code{window.size}. Since
#' the function won't perform calculations in that case, it might be faster, but if this behavior is not desired,
#' consider reinterpolating the time series (see \code{\link{reinterpolate}}) or increasing the window size.
#'
#' @references
#'
#' Cuturi, M. (2011). Fast global alignment kernels. In \emph{Proceedings of the 28th international conference on
#' machine learning (ICML-11)} (pp. 929-936).
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
#' dtwclust(CharTraj, k = 20L,
#'          distance = "gak", centroid = "shape",
#'          sigma = sigma,
#'          control = list(trace = TRUE,
#'                         window.size = 18L))
#' }
#'
GAK <- function(x, y, ..., sigma = NULL, window.size = NULL, normalize = TRUE) {
    x <- cbind(x)
    y <- cbind(y)

    consistency_check(x, "ts")
    consistency_check(y, "ts")

    ## check dimension consistency
    is_multivariate(list(x,y))

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
        window.size <- consistency_check(window.size, "window")

    logGAK <- .Call("logGAK", x, y,
                    NROW(x), NROW(y), NCOL(x),
                    sigma, window.size,
                    PACKAGE = "dtwclust")

    if (normalize) {
        gak_x <- .Call("logGAK", x, x, NROW(x), NROW(x), NCOL(x), sigma, window.size, PACKAGE = "dtwclust")
        gak_y <- .Call("logGAK", y, y, NROW(y), NROW(y), NCOL(y), sigma, window.size, PACKAGE = "dtwclust")

        logGAK <- 1 - exp(logGAK - 0.5 * (gak_x + gak_y))
    }

    attr(logGAK, "sigma") <- sigma

    logGAK
}

GAK_proxy <- function(x, y = NULL, ..., sigma = NULL, normalize = TRUE, pairwise = FALSE, symmetric = FALSE) {
    if (!normalize)
        warning("The proxy version of GAK is always normalized.")

    x <- consistency_check(x, "tsmat")
    consistency_check(x, "vltslist")

    dots <- list(...)

    ## normalization will be done manually to avoid multiple calculations of gak_x and gak_y
    dots$normalize <- FALSE

    if (is.null(y)) {
        y <- x
        symmetric <- TRUE

    } else {
        y <- consistency_check(y, "tsmat")
        consistency_check(y, "vltslist")

        if (symmetric && length(x) != length(y))
            stop("GAK: Dimension inconsistency when calculating cross-distance matrix within proxy.")
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

    gak_x <- sapply(x, function(x) do.call("GAK", enlist(x = x, y = x, dots = dots)))
    gak_y <- sapply(y, function(y) do.call("GAK", enlist(x = y, y = y, dots = dots)))

    retclass <- "crossdist"

    ## Register doSEQ if necessary
    check_parallel()

    ## Calculate distance matrix
    if (pairwise) {
        X <- split_parallel(x)
        Y <- split_parallel(y)

        validate_pairwise(X, Y)

        D <- foreach(x = X, y = Y,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %dopar% {
                         mapply(x, y, FUN = function(x, y) {
                             do.call("GAK",
                                     enlist(x = x,
                                            y = y,
                                            dots = dots))
                         })
                     }

        D <- 1 - exp(D - 0.5 * (gak_x + gak_y))

        names(D) <- NULL
        retclass <- "pairdist"

    } else if (symmetric) {
        pairs <- call_pairs(length(x), lower = FALSE)
        pairs <- split_parallel(pairs, 1L)

        d <- foreach(pairs = pairs,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %dopar% {
                         mapply(x[pairs[ , 1L]], x[pairs[ , 2L]],
                                SIMPLIFY = TRUE,
                                FUN = function(xx, yy) {
                                    do.call("GAK",
                                            enlist(x = xx,
                                                   y = yy,
                                                   dots = dots))
                                })
                     }

        rm("pairs")

        D <- matrix(0, nrow = length(x), ncol = length(x))
        D[upper.tri(D)] <- d
        D <- t(D)
        D[upper.tri(D)] <- d

        D <- 1 - exp(D - outer(gak_x, gak_y, function(x, y) { (x + y) / 2 }))
        diag(D) <- 0

        attr(D, "dimnames") <- list(names(x), names(x))

    } else {
        X <- split_parallel(x)

        D <- foreach(x = X,
                     .combine = rbind,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %dopar% {
                         ret <- lapply(x, y = y, FUN = function(x, y) {
                             sapply(y, x = x, FUN = function(y, x) {
                                 do.call("GAK",
                                         enlist(x = x,
                                                y = y,
                                                dots = dots))
                             })
                         })

                         do.call(rbind, ret)
                     }

        D <- 1 - exp(D - outer(gak_x, gak_y, function(x, y) { (x + y) / 2 }))
    }

    class(D) <- retclass
    attr(D, "method") <- "GAK"
    attr(D, "sigma") <- sigma

    D
}
