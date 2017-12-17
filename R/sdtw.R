#' Soft-DTW distance
#'
#' Soft-DTW distance measure as proposed in Cuturi and Blondel (2017).
#'
#' @export
#'
#' @param x,y Time series. Multivariate series must have time spanning the rows and variables
#'   spanning the columns.
#' @param gamma Positive regularization parameter, with lower values resulting in less smoothing.
#' @param ... Currently ignored.
#' @param cm Optionally, a matrix to use for the calculations. It should have `NROW(x)+1` rows and
#'   `NROW(y)+1` columns. Used internally for memory optimization. If provided, it **will** be
#'   modified *in place* by `C` code, except in the parallel version in [proxy::dist()] which
#'   ignores it for thread-safe reasons.
#' @template error-check
#'
#' @details
#'
#' Unlike other distances, soft-DTW can return negative values, and `sdtw(x, x)` is not always equal
#' to zero. Like DTW, soft-DTW does not follow the triangle inequality, but it is always symmetric.
#'
#' @return The Soft DTW distance.
#'
#' @template proxy
#' @template symmetric
#'
#' @references
#'
#' Cuturi, M., & Blondel, M. (2017). Soft-DTW: a Differentiable Loss Function for Time-Series. arXiv
#' preprint arXiv:1703.01541.
#'
sdtw <- function(x, y, gamma = 0.01, ..., cm = NULL, error.check = TRUE)
{
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }
    if (gamma <= 0) stop("The gamma paramter must be positive")
    mv <- is_multivariate(list(x,y)) # dimension consistency checked here

    if (is.null(cm))
        cm <- matrix(0, NROW(x) + 1L, NROW(y) + 1L)
    else if (!is.matrix(cm) || nrow(cm) < (NROW(x) + 1L) || ncol(cm) < (NROW(y) + 1L))
        stop("sdtw: Dimension inconsistency in 'cm'")
    else if (storage.mode(cm) != "double")
        stop("sdtw: If provided, 'cm' must have 'double' storage mode.")

    # return
    .Call(C_soft_dtw, x, y, gamma, cm, NULL, mv, PACKAGE = "dtwclust")
}

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

#' @importFrom bigmemory attach.big.matrix
#' @importFrom bigmemory describe
#' @importFrom bigmemory is.big.matrix
#'
sdtw_proxy <- function(x, y = NULL, ..., cm = NULL, error.check = TRUE, pairwise = FALSE) {
    x <- tslist(x)
    if (error.check) check_consistency(x, "vltslist")

    dots <- list(...)
    retclass <- "crossdist"

    if (is.null(y)) {
        y <- x
        symmetric <- TRUE

    } else {
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
        symmetric <- FALSE
    }

    # pre-allocate cm
    if(is.null(cm)) cm <- matrix(0, max(sapply(x, NROW)) + 1L, max(sapply(y, NROW)) + 1L)
    pairwise <- isTRUE(pairwise)
    dim_out <- c(length(x), length(y))
    dim_names <- list(names(x), names(y))
    D <- allocate_distmat(length(x), length(y), pairwise, symmetric) # utils.R

    # Wrap as needed for foreach
    if (pairwise) {
        x <- split_parallel(x)
        y <- split_parallel(y)
        validate_pairwise(x, y)
        endpoints <- attr(x, "endpoints")

    } else if (symmetric) {
        endpoints <- symmetric_loop_endpoints(length(x)) # utils.R
        diagonal <- sdtw_proxy(x, ..., cm = cm, error.check = FALSE, pairwise = TRUE)
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- x

    } else {
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- split_parallel(y)
        endpoints <- attr(y, "endpoints")
    }

    if (bigmemory::is.big.matrix(D)) {
        D_desc <- bigmemory::describe(D)
        noexport <- "D"
        packages <- c("dtwclust", "bigmemory")

    } else {
        D_desc <- NULL
        noexport <- ""
        packages <- c("dtwclust")
    }

    # Calculate distance matrix
    foreach(x = x, y = y, endpoints = endpoints,
            .combine = c,
            .multicombine = TRUE,
            .packages = packages,
            .export = c("sdtw_loop", "enlist"),
            .noexport = noexport) %op% {
                bigmat <- !is.null(D_desc)
                d <- if (bigmat) bigmemory::attach.big.matrix(D_desc)@address else D
                do.call(sdtw_loop,
                        enlist(d = d,
                               x = x,
                               y = y,
                               symmetric = symmetric,
                               pairwise = pairwise,
                               endpoints = endpoints,
                               bigmat = bigmat,
                               cm = cm,
                               dots = dots),
                        TRUE)
            }

    D <- D[,]
    if (pairwise) {
        class(D) <- "pairdist"

    } else {
        if (is.null(dim(D))) dim(D) <- dim_out
        dimnames(D) <- dim_names
        if (symmetric) D[cbind(1L:dim_out[1L], 1L:dim_out[2L])] <- diagonal
        class(D) <- "crossdist"
    }

    attr(D, "method") <- "SDTW"
    # return
    D
}

# ==================================================================================================
# Wrapper for C++
# ==================================================================================================

sdtw_loop <- function(d, x, y, symmetric, pairwise, endpoints, bigmat, ...,
                      gamma = 0.01, cm = NULL)
{
    mv <- is_multivariate(c(x, y))
    if (gamma <= 0) stop("The gamma paramter must be positive")

    nr <- max(sapply(x, NROW)) + 1L
    nc <- max(sapply(y, NROW)) + 1L
    if (!is.matrix(cm) || nrow(cm) < nr || ncol(cm) < nc)
        stop("sdtw: Dimension inconsistency in 'cm'")
    if (storage.mode(cm) != "double")
        stop("sdtw: If provided, 'cm' must have 'double' storage mode.")

    distargs <- list(gamma = gamma, cm = cm)

    .Call(C_sdtw_loop,
          d, x, y, distargs, symmetric, pairwise, bigmat, mv, endpoints,
          PACKAGE = "dtwclust")
}
