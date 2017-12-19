sdtw_cent_nloptr <- function(centroid, series, gamma, weights, mv, dim0, cm, dm, em)
{
    if (mv && is.null(dim(centroid))) dim(centroid) <- dim0
    .Call(C_sdtw_cent, series, centroid, gamma, weights, mv, cm, dm, em, PACKAGE = "dtwclust")
}

#' Centroid calculation based on soft-DTW
#'
#' Soft-DTW centroid function as proposed in Cuturi and Blondel (2017).
#'
#' @export
#' @importFrom nloptr nloptr
#'
#' @param series A matrix or data frame where each row is a time series, or a list where each
#'   element is a time series. Multivariate series should be provided as a list of matrices where
#'   time spans the rows and the variables span the columns of each matrix.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of
#'   `series` if `NULL`. For multivariate series, this should be a matrix with the same
#'   characteristics as the matrices in `series`.
#' @param gamma Positive regularization parameter, with lower values resulting in less smoothing.
#' @param weights A vector of weights for each element of `series`.
#' @param ... Further arguments for [nloptr::nloptr()] (except `opts`).
#' @param opts List of options to pass to [nloptr::nloptr()].
#' @template error-check
#' @param cm,dm,em Optional helper matrices for the calculations. Used internally for memory
#'   optimization. If provided, they **will** be modified *in place* by `C` code. See details for
#'   dimensioning information.
#'
#' @details
#'
#' The helper matrices are allocated in the following way by default:
#'
#' ```
#' num_rows <- NROW(centroid)
#' num_cols <- max(sapply(series, NROW))
#' cm <- matrix(0, num_rows + 2L, num_cols + 2L)
#' dm <- matrix(0, num_rows + 1L, num_cols + 1L)
#' em <- matrix(0, 2L, num_cols + 2L)
#' ```
#'
#' Note that you can trace the optimization by specifing `print_level > 0` in `opts`.
#'
#' @return The resulting centroid, with attribute `nloptr_results` specifying the optimization
#' results (except for `solution`, which is the returned centroid).
#'
#' @references
#'
#' Cuturi, M., & Blondel, M. (2017). Soft-DTW: a Differentiable Loss Function for Time-Series. arXiv
#' preprint arXiv:1703.01541.
#'
sdtw_cent <- function(series, centroid = NULL, gamma = 0.01, weights = rep(1, length(series)), ...,
                      opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 20L),
                      error.check = TRUE, cm = NULL, dm = NULL, em = NULL)
{
    series <- tslist(series)
    if (is.null(centroid)) centroid <- series[[sample(length(series), 1L)]] # Random choice
    if (gamma <= 0) stop("The gamma paramter must be positive")
    mv <- is_multivariate(c(series, list(centroid)))
    if (length(weights) != length(series))
        stop("The 'weights' vector must have the same length as 'series'")
    if (error.check) {
        check_consistency(series, "vltslist")
        check_consistency(centroid, "ts")
    }

    # for helper matrices
    num_rows <- NROW(centroid)
    num_cols <- max(sapply(series, NROW))

    if (is.null(cm))
        cm <- matrix(0, num_rows + 2L, num_cols + 2L)
    else if (!is.matrix(cm) || nrow(cm) < (num_rows + 2L) || ncol(cm) < (num_cols + 2L))
        stop("sdtw_cent: Dimension inconsistency in 'cm'")
    else if (storage.mode(cm) != "double")
        stop("sdtw_cent: If provided, 'cm' must have 'double' storage mode.")

    if (is.null(dm))
        dm <- matrix(0, num_rows + 1L, num_cols + 1L)
    else if (!is.matrix(dm) || nrow(dm) < (num_rows + 1L) || ncol(dm) < (num_cols + 1L))
        stop("sdtw_cent: Dimension inconsistency in 'dm'")
    else if (storage.mode(dm) != "double")
        stop("sdtw_cent: If provided, 'dm' must have 'double' storage mode.")

    if (is.null(em))
        em <- matrix(0, 2L, num_cols + 2L)
    else if (!is.matrix(em) || nrow(em) < 2L || ncol(em) < (num_cols + 2L))
        stop("sdtw_cent: Dimension inconsistency in 'em'")
    else if (storage.mode(em) != "double")
        stop("sdtw_cent: If provided, 'em' must have 'double' storage mode.")

    dim0 <- dim(centroid)
    if (mv) nm0 <- dimnames(centroid)
    opt <- nloptr::nloptr(centroid, sdtw_cent_nloptr, opts = opts,
                          series = series, gamma = gamma, weights = weights, mv = mv, dim0 = dim0,
                          cm = cm, dm = dm, em = em)

    cent_out <- opt$solution
    opt$solution <- NULL
    if (mv) {
        dim(cent_out) <- dim0
        dimnames(cent_out) <- nm0
    }
    attr(cent_out, "nloptr_results") <- opt
    # return
    cent_out
}
