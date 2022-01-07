sdtw_cent_nloptr <- function(centroid, series, gamma, weights, mv, dim0) {
    if (mv && is.null(dim(centroid))) dim(centroid) <- dim0
    num_threads <- get_nthreads()
    .Call(C_sdtw_cent, series, centroid, gamma, weights, mv, num_threads, PACKAGE = "dtwclust")
}

sdtw_cent_stats <- function(centroid, series, gamma, weights, mv, dim0, .shared_env_) {
    objective_and_gradient <- sdtw_cent_nloptr(
        centroid,
        series,
        gamma,
        weights,
        mv,
        dim0
    )

    .shared_env_$gradient <- objective_and_gradient$gradient
    objective_and_gradient$objective
}

sdtw_cent_stats_gr <- function(ignored, .shared_env_, ...) {
    .shared_env_$gradient
}

#' Centroid calculation based on soft-DTW
#'
#' Soft-DTW centroid function as proposed in Cuturi and Blondel (2017).
#'
#' @export
#' @importFrom stats optim
#'
#' @param series A matrix or data frame where each row is a time series, or a list where each
#'   element is a time series. Multivariate series should be provided as a list of matrices where
#'   time spans the rows and the variables span the columns of each matrix.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of
#'   `series` if `NULL`. For multivariate series, this should be a matrix with the same
#'   characteristics as the matrices in `series`.
#' @param gamma Positive regularization parameter, with lower values resulting in less smoothing.
#' @param weights A vector of weights for each element of `series`.
#' @param ... Further arguments for the optimization backend (except `opts` for `nloptr`, `control`
#'   for `optim`, and `...` for both).
#' @param opts List of options to pass to `nloptr` or [stats::optim()]'s `control`. The defaults in
#'   the function's formals are for `nloptr`, but the value will be adjusted for `optim` if needed.
#' @template error-check
#'
#' @details
#'
#' This function can delegate the optimization to the \pkg{nloptr} package. For that to happen, you
#' must load it with either [base::library()] or [base::loadNamespace()]. If the aforementioned is
#' not fulfilled, the function will delegate to [stats::optim()].
#'
#' @return The resulting centroid, with the optimization results as attributes (except for the
#'   returned centroid).
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
#' @references
#'
#' Cuturi, M., & Blondel, M. (2017). Soft-DTW: a Differentiable Loss Function for Time-Series. arXiv
#' preprint arXiv:1703.01541.
#'
sdtw_cent <- function(series, centroid = NULL, gamma = 0.01, weights = rep(1, length(series)), ...,
                      error.check = TRUE,
                      opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 20L))
{
    series <- tslist(series)

    if (is.null(centroid)) centroid <- series[[sample(length(series), 1L)]] # Random choice
    if (gamma <= 0) stop("The gamma paramter must be positive")

    mv <- is_multivariate(c(series, list(centroid)))

    if (length(weights) != length(series)) {
        stop("The 'weights' vector must have the same length as 'series'")
    }
    if (error.check) {
        check_consistency(series, "vltslist")
        check_consistency(centroid, "ts")
    }

    gamma <- as.numeric(gamma)[1L]
    weights <- as.numeric(weights)

    dim0 <- dim(centroid)
    if (mv) nm0 <- dimnames(centroid)

    dots <- list(...)

    if (isNamespaceLoaded("nloptr")) { # nocov start
        backend <- "nloptr"
        nloptr <- get("nloptr", asNamespace("nloptr"), mode = "function")
        dots <- dots[intersect(names(dots), setdiff(names(formals(nloptr)), "..."))]

        opt <- quoted_call(
            nloptr,
            x0 = centroid,
            eval_f = sdtw_cent_nloptr,
            opts = opts,
            series = series,
            gamma = gamma,
            weights = weights,
            mv = mv,
            dim0 = dim0,
            dots = dots
        )

        cent_out <- opt$solution
    } # nocov end
    else {
        backend <- "stats"
        dots <- dots[intersect(names(dots), setdiff(names(formals(stats::optim)), "..."))]

        if (missing(opts)) {
            opts <- list(maxit = 20L)
        }

        opt <- quoted_call(
            stats::optim,
            par = centroid,
            fn = sdtw_cent_stats,
            gr = sdtw_cent_stats_gr,
            method = "L-BFGS-B",
            control = opts,
            series = series,
            gamma = gamma,
            weights = weights,
            mv = mv,
            dim0 = dim0,
            .shared_env_ = new.env(),
            dots = dots
        )

        cent_out <- opt$par
    }

    if (mv) {
        dim(cent_out) <- dim0
        dimnames(cent_out) <- nm0
    }

    if (getOption("dtwclust_sdtw_cent_return_attrs", TRUE)) { # nocov start
        opt$par <- opt$call <- opt$solution <- NULL
        attr(cent_out, paste0(backend, "_results")) <- opt
    } # nocov end

    # return
    cent_out
}
