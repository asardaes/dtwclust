#' DTW Barycenter Averaging
#'
#' A global averaging method for time series under DTW (Petitjean, Ketterlin and Gancarski 2011).
#'
#' @export
#'
#' @param X A matrix or data frame where each row is a time series, or a list where each element is
#'   a time series. Multivariate series should be provided as a list of matrices where time spans
#'   the rows and the variables span the columns of each matrix.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of `X`
#'   if `NULL`. For multivariate series, this should be a matrix with the same characteristics as
#'   the matrices in `X`.
#' @param ... Further arguments for [dtw_basic()]. However, the following are already pre-
#'   specified: `window.size`, `norm` (passed along), and `backtrack`.
#' @param window.size Window constraint for the DTW calculations. `NULL` means no constraint. A
#'   slanted band is used by default.
#' @param norm Norm for the local cost matrix of DTW. Either "L1" for Manhattan distance or "L2" for
#'   Euclidean distance.
#' @param max.iter Maximum number of iterations allowed.
#' @param delta At iteration `i`, if `all(abs(centroid_{i}` `-` `centroid_{i-1})` `< delta)`,
#'   convergence is assumed.
#' @template error-check
#' @param trace If `TRUE`, the current iteration is printed to output.
#' @param gcm Optional matrix to pass to [dtw_basic()] (for the case when `backtrack = TRUE`). To
#'   define the matrix size, it should be assumed that `x` is the *longest* series in `X`, and `y`
#'   is the `centroid` if provided or `x` otherwise.
#' @param mv.ver Multivariate version to use. See below.
#'
#' @details
#'
#' This function tries to find the optimum average series between a group of time series in DTW
#' space. Refer to the cited article for specific details on the algorithm.
#'
#' If a given series reference is provided in `centroid`, the algorithm should always converge to
#' the same result provided the elements of `X` keep the same values, although their order may
#' change.
#'
#' @template window
#'
#' @return The average time series.
#'
#' @section Multivariate series:
#'
#'   There are currently 2 versions of DBA implemented for multivariate series:
#'
#'   - If `mv.ver = "by-variable"`, then each variable of each series in `X` and `centroid` are
#'     extracted, and the univariate version of the algorithm is applied to each set of variables,
#'     binding the results by column. Therefore, the DTW backtracking is different for each
#'     variable.
#'   - If `mv.ver = "by-series"`, then all variables are considered at the same time, so the DTW
#'     backtracking is computed based on each multivariate series as a whole. This version was
#'     implemented in version 4.0.0 of \pkg{dtwclust}, and it might be faster.
#'
#' @note
#'
#' The indices of the DTW alignment are obtained by calling [dtw_basic()] with `backtrack = TRUE`.
#'
#' @references
#'
#' Petitjean F, Ketterlin A and Gancarski P (2011). ``A global averaging method for dynamic time
#' warping, with applications to clustering.'' *Pattern Recognition*, **44**(3), pp. 678 - 693. ISSN
#' 0031-3203, \url{http://dx.doi.org/10.1016/j.patcog.2010.09.013},
#' \url{http://www.sciencedirect.com/science/article/pii/S003132031000453X}.
#'
#' @examples
#'
#' # Sample data
#' data(uciCT)
#'
#' # Obtain an average for the first 5 time series
#' dtw_avg <- DBA(CharTraj[1:5], CharTraj[[1]], trace = TRUE)
#'
#' # Plot
#' matplot(do.call(cbind, CharTraj[1:5]), type = "l")
#' points(dtw_avg)
#'
#' # Change the provided order
#' dtw_avg2 <- DBA(CharTraj[5:1], CharTraj[[1]], trace = TRUE)
#'
#' # Same result?
#' all(dtw_avg2 == dtw_avg2)
#'
DBA <- function(X, centroid = NULL, ...,
                window.size = NULL, norm = "L1",
                max.iter = 20L, delta = 1e-3,
                error.check = TRUE, trace = FALSE,
                gcm = NULL, mv.ver = "by-variable")
{
    X <- any2list(X)
    mv.ver <- match.arg(mv.ver, c("by-variable", "by-series"))
    mv.ver <- switch(mv.ver, "by-variable" = 1L, "by-series" = 2L)

    if (is.null(centroid)) centroid <- X[[sample(length(X), 1L)]] # Random choice

    if (error.check) {
        check_consistency(X, "vltslist")
        check_consistency(centroid, "ts")
    }

    if (is.null(window.size))
        window.size <- -1L
    else
        window.size <- check_consistency(window.size, "window")

    if (max.iter < 1L)
        stop("Maximum iterations must be positive.")
    else
        max.iter <- as.integer(max.iter)[1L]

    dots <- list(...)
    step.pattern <- dots$step.pattern

    if (is.null(step.pattern) || identical(step.pattern, symmetric2))
        step.pattern <- 2
    else if (identical(step.pattern, symmetric1))
        step.pattern <- 1
    else
        stop("step.pattern must be either symmetric1 or symmetric2 (without quotes)")

    if (length(delta) > 1L) delta <- delta[1L]
    trace <- isTRUE(trace)
    norm <- match.arg(norm, c("L1", "L2"))
    norm <- switch(norm, "L1" = 1, "L2" = 2)
    L <- max(sapply(X, NROW)) + 1L ## maximum length of considered series + 1L
    mv <- is_multivariate(c(X, list(centroid)))
    nr <- NROW(centroid) + 1L

    ## pre-allocate cost matrices
    if (is.null(gcm))
        gcm <- matrix(0, L, nr)
    else if (!is.matrix(gcm) || nrow(gcm) < L || ncol(gcm) < nr)
        stop("DBA: Dimension inconsistency in 'gcm'")
    else if (storage.mode(gcm) != "double")
        stop("DBA: If provided, 'gcm' must have 'double' storage mode.")

    ## All parameters for dtw_basic()
    dots <- list(window.size = window.size,
                 norm = norm,
                 gcm = gcm,
                 step.pattern = step.pattern,
                 backtrack = TRUE)

    ## C++ code
    new_cent <- .Call(C_dba, X, centroid, max.iter, delta, trace, mv, mv.ver, dots, PACKAGE = "dtwclust")
    if (mv) dimnames(new_cent) <- dimnames(centroid)
    new_cent
}
