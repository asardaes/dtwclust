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
#' @template error-check
#'
#' @details
#'
#' Unlike other distances, soft-DTW can return negative values, and `sdtw(x, x)` is not always equal
#' to zero. Like DTW, soft-DTW does not fulfill the triangle inequality, but it is always symmetric.
#'
#' @return The Soft DTW distance.
#'
#' @template rcpp-parallel
#' @template proxy
#' @template symmetric
#'
#' @references
#'
#' Cuturi, M., & Blondel, M. (2017). Soft-DTW: a Differentiable Loss Function for Time-Series. arXiv
#' preprint arXiv:1703.01541.
#'
sdtw <- function(x, y, gamma = 0.01, ..., error.check = TRUE)
{
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
    if (error.check) check_consistency(x, "vltslist")
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
    if (!pairwise && symmetric)
        diagonal <- sdtw_proxy(x, gamma = gamma, error.check = FALSE, pairwise = TRUE)
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
    else {
        dimnames(D) <- dim_names
        if (symmetric) diag(D) <- diagonal
        class(D) <- "crossdist"
    }
    attr(D, "method") <- "SDTW"
    # return
    D
}
