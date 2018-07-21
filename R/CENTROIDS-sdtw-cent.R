sdtw_cent_nloptr <- function(centroid, series, gamma, weights, mv, dim0)
{
    if (mv && is.null(dim(centroid))) dim(centroid) <- dim0
    num_threads <- get_nthreads()
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
#'
#' @details
#'
#' Note that you can trace the optimization by specifying `print_level > 0` in `opts`.
#'
#' @template rcpp-parallel
#'
#' @section Parallel Computing:
#'
#'   For unknown reasons, this function has returned different results (in the order of 1e-6) when
#'   using multi-threading in x64 Windows installations in comparison to other environments (using
#'   nloptr v1.0.4). Consider limiting the number of threads if you run into reproducibility
#'   problems.
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
                      error.check = TRUE)
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
        dim0 = dim0
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
