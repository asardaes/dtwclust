#' TADPole clustering
#'
#' Time-series Anytime Density Peaks Clustering as proposed by Begum et al. (2015).
#'
#' @export
#'
#' @param data A matrix or data frame where each row is a time series, or a list where each element
#'   is a time series. Multivariate series are **not** supported.
#' @param window.size Window size constraint for DTW (Sakoe-Chiba). See details.
#' @param k The number of desired clusters. Can be a vector with several values.
#' @param dc The cutoff distance(s). May be a vector with several values.
#' @template error-check
#' @param lb Which lower bound to use, "lbk" for [lb_keogh()] or "lbi" for [lb_improved()].
#' @param trace Logical flag. If `TRUE`, more output regarding the progress is printed to screen.
#'
#' @details
#'
#' This function can be called either directly or through [dtwclust()] and [tsclust()].
#'
#' TADPole clustering adopts a relatively new clustering framework and adapts it to time series
#' clustering with DTW. See the cited article for the details of the algorithm.
#'
#' Because of the way the algorithm works, it can be considered a kind of Partitioning Around
#' Medoids (PAM). This means that the cluster centroids are always elements of the data. However,
#' this algorithm is deterministic, depending on the value of `dc`.
#'
#' The algorithm first uses the DTW's upper and lower bounds (Euclidean and LB_Keogh respectively)
#' to find series with many close neighbors (in DTW space). Anything below the cutoff distance
#' (`dc`) is considered a neighbor. Aided with this information, the algorithm then tries to prune
#' as many DTW calculations as possible in order to accelerate the clustering procedure. The series
#' that lie in dense areas (i.e. that have lots of neighbors) are taken as cluster centroids.
#'
#' The algorithm relies on the DTW bounds, which are only defined for univariate time series of
#' equal length.
#'
#' Parallelization is supported, but given the internal optimizations, it may only be useful if
#' multiple `dc` values are specified in the same call.
#'
#' @template window
#'
#' @return
#'
#' A list with:
#'
#' - `cl`: Cluster indices.
#' - `centroids`: Indices of the centroids.
#' - `distCalcPercentage`: Percentage of distance calculations that were actually performed.
#'
#' For multiple `k`/`dc` values, a list of lists is returned, each internal list having the
#' aforementioned elements.
#'
#' @references
#'
#' Begum N, Ulanova L, Wang J and Keogh E (2015). ``Accelerating Dynamic Time Warping Clustering
#' with a Novel Admissible Pruning Strategy.'' In *Conference on Knowledge Discovery and Data
#' Mining*, series KDD '15. ISBN 978-1-4503-3664-2/15/08, \url{
#' http://dx.doi.org/10.1145/2783258.2783286}.
#'
TADPole <- function(data, k = 2L, dc, window.size, error.check = TRUE, lb = "lbk", trace = FALSE) {
    if (missing(window.size)) stop("Please provide a positive window size")
    if (missing(dc)) stop("Please provide the 'dc' parameter")
    if (any(dc < 0)) stop("The cutoff distance 'dc' must be positive")

    x <- any2list(data)
    n <- length(x)

    if (n < 2L) stop("data should have more than one time series")
    if (any(k > n)) stop("Number of clusters should be less than the number of time series")
    lb <- match.arg(lb, c("lbk", "lbi"))

    if (trace) cat("\tComputing lower and upper bound matrices\n\n")

    ## Calculate matrices with bounds (error check in lbk/lbi)
    LBM <- proxy::dist(x, x,
                       method = lb,
                       window.size = window.size,
                       force.symmetry = TRUE,
                       norm = "L2",
                       error.check = error.check)

    ## NOTE: Euclidean is only valid as upper bound if 'symmetric1' step pattern is used
    UBM <- proxy::dist(x, x, method = "L2")

    len <- max(lengths(data)) + 1L
    dtw_args <- list(window.size = window.size,
                     norm = 2,
                     step.pattern = 1,
                     backtrack = FALSE,
                     gcm = matrix(0, 2L, len))

    RET <- foreach(
        dc = dc, .combine = c, .multicombine = TRUE,
        .packages = "dtwclust", .export = "call_tadpole") %op%
        {
            ret <- lapply(k, function(dummy) { list() }) ## modified in place in C++
            call_tadpole(x, k, dc, dtw_args, LBM, UBM, trace, ret) ## foreach can't see C objects
            ## return
            ret
        }

    ## Return
    if (length(RET) == 1L) RET[[1L]] else RET
}

## a parallel foreach loop will not see the C objects, so I need to pass a function from the
## dtwclust namespace (via the .export parameter)
call_tadpole <- function(x, k, dc, dtw_args, LBM, UBM, trace, ret) {
    .Call(C_tadpole, x, k, dc, dtw_args, LBM, UBM, trace, ret, PACKAGE = "dtwclust")
}
