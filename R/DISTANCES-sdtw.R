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
#' @param error.check `r roxygen_error_check_param()`
#'
#' @details
#'
#' Unlike other distances, soft-DTW can return negative values, and `sdtw(x, x)` is not always equal
#' to zero. Like DTW, soft-DTW does not fulfill the triangle inequality, but it is always symmetric.
#'
#' @return The Soft DTW distance.
#'
#' @section `r roxygen_proxy_section()`
#'
#' `r roxygen_proxy_symmetric()`
#'
#' Note that, due to the fact that this distance is not always zero when a series is compared
#' against itself, this optimization is likely problematic for soft-DTW, as the `dist` object will
#' be handled by many functions as if it had only zeroes in the diagonal. An exception is
#' [tsclust()] when using partitional clustering with PAM centroids---actual diagonal values will
#' be calculated and considered internally in that case.
#'
#' @references
#'
#' Cuturi, M., & Blondel, M. (2017). Soft-DTW: a Differentiable Loss Function for Time-Series. arXiv
#' preprint arXiv:1703.01541.
#'
sdtw <- function(x, y, gamma = 0.01, ..., error.check = TRUE) {
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }
    if (gamma <= 0) stop("The gamma paramter must be positive")
    mv <- is_multivariate(list(x,y)) # dimension consistency checked here
    cm <- matrix(0, NROW(x) + 2L, NROW(y) + 2L)
    # return
    .Call(C_soft_dtw, x, y, gamma, cm, mv, PACKAGE = "dtwclust")
}

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

sdtw_proxy <- function(x, y = NULL, gamma = 0.01, ..., error.check = TRUE, pairwise = FALSE) {
    x <- tslist(x)
    if (error.check) {
        check_consistency(x, "vltslist")
    }

    if (is.null(y)) {
        y <- x
        symmetric <- TRUE
    }
    else {
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
        symmetric <- FALSE
    }

    fill_type <- mat_type <- dim_names <- NULL # avoid warning about undefined globals
    eval(prepare_expr) # UTILS-expressions.R

    # adjust parameters for this distance
    if (gamma <= 0) stop("The 'gamma' parameter must be positive")
    mv <- is_multivariate(c(x, y))

    # calculate distance matrix
    distance <- "SDTW" # read in C++, can't be temporary!
    distargs <- list(
        gamma = gamma
    )
    num_threads <- get_nthreads()
    .Call(C_distmat_loop,
          D, x, y, distance, distargs, fill_type, mat_type, num_threads,
          PACKAGE = "dtwclust")

    if (pairwise) {
        dim(D) <- NULL
        class(D) <- "pairdist"
    }
    else if (symmetric) {
        warning("The distance between a series and itself is not always 0 with soft-DTW,",
                " and this will be hidden in the 'dist' object that only includes lower triangular values.")
        dim(D) <- NULL
        class(D) <- "dist"
        attr(D, "Size") <- length(x)
        attr(D, "Labels") <- names(x)
    }
    else {
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }

    attr(D, "method") <- "SDTW"
    D
}

sdtw_wrapper <- function(x, gamma = 0.01, ...) {
    x <- tslist(x)
    check_consistency(x, "vltslist")

    D <- allocate_distmat(length(x), length(x), FALSE, TRUE, TRUE) # UTILS-utils.R
    fill_type <- "LOWER_TRIANGULAR_DIAGONAL"
    mat_type <- "R_MATRIX"

    # adjust parameters for this distance
    if (gamma <= 0) stop("The 'gamma' parameter must be positive")
    mv <- is_multivariate(x)

    # calculate distance matrix
    distance <- "SDTW" # read in C++, can't be temporary!
    distargs <- list(
        gamma = gamma
    )
    num_threads <- get_nthreads()

    .Call(C_distmat_loop,
          D, x, x, distance, distargs, fill_type, mat_type, num_threads,
          PACKAGE = "dtwclust")

    dim(D) <- NULL
    class(D) <- "dist"
    attr(D, "Size") <- length(x)
    attr(D, "Labels") <- names(x)
    attr(D, "Diag") <- TRUE
    attr(D, "method") <- "SDTW"

    D
}
