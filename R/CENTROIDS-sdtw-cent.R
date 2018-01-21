sdtw_cent_nloptr <- function(centroid, series, gamma, weights, mv, dim0, num_threads)
{
    if (mv && is.null(dim(centroid))) dim(centroid) <- dim0
    .Call(C_sdtw_cent, series, centroid, gamma, weights, mv, num_threads, PACKAGE = "dtwclust")
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
#' @param ... Further arguments for [nloptr::nloptr()] (except `opts` and `...`).
#' @param opts List of options to pass to [nloptr::nloptr()].
#' @template error-check
#' @param num_threads How many threads to use for multi-threading. See Parallel Computing section
#'   below.
#'
#' @details
#'
#' Note that you can trace the optimization by specifing `print_level > 0` in `opts`.
#'
#' @template rcpp-parallel
#'
#' @section Parallel Computing:
#'
#'   In contrast to other \pkg{dtwclust} functions, this function has a parameter to specify how
#'   many threads it should try to use. This is because the procedure is more sensitive to floating
#'   point inaccuracies, and dividing the work onto different threads could change the results (at
#'   least in the order of 1e-8 according to my limited experiments). This could be an issue during
#'   clusterings where a lot of calls to the centroid function are made and the errors propagate.
#'   Note that you should still set the maximum amount of available threads with
#'   [RcppParallel::setThreadOptions()].
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
                      error.check = TRUE, num_threads = 1L)
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
    gamma <- as.numeric(gamma)[1L]
    weights <- as.numeric(weights)
    num_threads <- as.integer(num_threads)[1L]

    dots <- list(...)
    dots <- dots[intersect(names(dots), setdiff(names(formals(nloptr::nloptr)), "..."))]
    dim0 <- dim(centroid)
    if (mv) nm0 <- dimnames(centroid)
    opt <- do.call(what = nloptr::nloptr, quote = TRUE, args = enlist(
        x0 = centroid,
        eval_f = sdtw_cent_nloptr,
        opts = opts,
        dots = dots,
        series = series,
        gamma = gamma,
        weights = weights,
        mv = mv,
        dim0 = dim0,
        num_threads = num_threads
    ))

    cent_out <- opt$solution
    opt$call <- opt$solution <- NULL
    if (mv) {
        dim(cent_out) <- dim0
        dimnames(cent_out) <- nm0
    }
    attr(cent_out, "nloptr_results") <- opt
    # return
    cent_out
}
