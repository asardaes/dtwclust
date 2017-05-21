#' DTW distance matrix guided by Lemire's improved lower bound
#'
#' Calculation of a distance matrix with the Dynamic Time Warping (DTW) distance guided by Lemire's
#' improved lower bound (LB_Improved).
#'
#' @export
#'
#' @param x,y A matrix or data frame where rows are time series, or a list of time series.
#' @param window.size Window size to use with the LB and DTW calculation. See details.
#' @param norm Pointwise distance. Either `"L1"` for Manhattan distance or `"L2"` for Euclidean.
#' @template error-check
#' @param pairwise Calculate pairwise distances?
#' @param dtw.func Which function to use for core DTW the calculations, either "dtw" or "dtw_basic".
#'   See [dtw::dtw()] and [dtw_basic()].
#' @param ... Further arguments for `dtw.func` or [lb_improved()].
#'
#' @details
#'
#' This function first calculates an initial estimate of a distance matrix between two sets of time
#' series using LB_Improved. Afterwards, it uses the estimate to calculate the corresponding true
#' DTW distance between *only* the nearest neighbors of each series in `x` found in `y`, and it
#' continues iteratively until no changes in the nearest neighbors occur.
#'
#' If only `x` is provided, the distance matrix is calculated between all its time series.
#'
#' This could be useful in case one is interested in only the nearest neighbor of one or more series
#' within a dataset.
#'
#' @template window
#'
#' @return The distance matrix with class `crossdist`.
#'
#' @template parallel
#'
#' @note
#'
#' This function uses a lower bound that is only defined for time series of equal length.
#'
#' A considerably large dataset is probably necessary before this is faster than using [dtw_basic()]
#' with [proxy::dist()].
#'
#' Nearest neighbors are found row-wise.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .'' *Pattern
#' Recognition*, **42**(9), pp. 2169 - 2180. ISSN 0031-3203,
#' \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030},
#' \url{http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
#'
#' @seealso
#'
#' [lb_keogh()], [lb_improved()]
#'
#' @examples
#'
#' # Load data
#' data(uciCT)
#'
#' # Reinterpolate to same length
#' data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
#'
#' # Calculate the DTW distance between a certain subset aided with the lower bound
#' system.time(d <- dtw_lb(data[1:5], data[6:50], window.size = 20))
#'
#' # Nearest neighbors
#' NN1 <- apply(d, 1L, which.min)
#'
#' # Calculate the DTW distances between all elements (slower)
#' system.time(d2 <- proxy::dist(data[1:5], data[6:50], method = "DTW",
#'                               window.type = "slantedband", window.size = 20))
#'
#' # Nearest neighbors
#' NN2 <- apply(d2, 1L, which.min)
#'
#' # Calculate the DTW distances between all elements using dtw_basic (actually faster, see notes)
#' system.time(d3 <- proxy::dist(data[1:5], data[6:50], method = "DTW_BASIC",
#'                               window.size = 20))
#'
#' # Nearest neighbors
#' NN3 <- apply(d3, 1L, which.min)
#'
#' # Same results?
#' all(NN1 == NN2)
#' all(NN1 == NN3)
#'
#' \dontrun{
#' #### Running DTW_LB with parallel support
#' # For such a small dataset, this is probably slower in parallel
#' require(doParallel)
#'
#' # Create parallel workers
#' cl <- makeCluster(detectCores())
#' invisible(clusterEvalQ(cl, library(dtwclust)))
#' registerDoParallel(cl)
#'
#' # Distance matrix
#' D <- dtw_lb(data[1:50], data[51:100], window.size = 20)
#'
#' # Stop parallel workers
#' stopCluster(cl)
#'
#' # Return to sequential computations
#' registerDoSEQ()
#'
#' # Nearest neighbors
#' NN <- apply(D, 1, which.min)
#' cbind(names(data[1:50]), names(data[51:100][NN]))
#' }
#'
dtw_lb <- function(x, y = NULL, window.size = NULL, norm = "L1",
                   error.check = TRUE, pairwise = FALSE,
                   dtw.func = "dtw_basic", ...)
{
    norm <- match.arg(norm, c("L1", "L2"))
    dtw.func <- match.arg(dtw.func, c("dtw", "dtw_basic"))

    if (dtw.func == "dtw")
        method <- if (norm == "L1") "DTW" else "DTW2"
    else
        method <- toupper(dtw.func)

    X <- any2list(x)

    if (is.null(y))
        Y <- X
    else
        Y <- any2list(y)

    if (is_multivariate(X) || is_multivariate(Y))
        stop("dtw_lb does not support multivariate series.")

    dots <- list(...)
    dots$dist.method <- norm
    dots$norm <- norm
    dots$window.size <- window.size
    dots$pairwise <- TRUE

    if (pairwise) {
        check_consistency(X, "tslist")
        check_consistency(Y, "tslist")

        if (is.null(window.size))
            dots$window.type <- "none"
        else
            dots$window.type <- "slantedband"

        X <- split_parallel(X)
        Y <- split_parallel(Y)
        validate_pairwise(X, Y)

        D <- foreach(X = X, Y = Y,
                     .combine = c,
                     .multicombine = TRUE,
                     .packages = "dtwclust",
                     .export = "enlist") %op% {
                         do.call(proxy::dist,
                                 enlist(x = X, y = Y,
                                        method = method,
                                        dots = dots))
                     }

        return(D)
    }

    window.size <- check_consistency(window.size, "window")
    dots$window.size <- window.size
    dots$window.type <- "slantedband"

    ## NOTE: I tried starting with LBK estimate, refining with LBI and then DTW but, overall,
    ## it was usually slower, almost the whole matrix had to be recomputed for LBI

    ## Initial estimate
    D <- proxy::dist(X, Y, method = "LBI",
                     window.size = window.size,
                     norm = norm,
                     error.check = error.check,
                     ...)

    D <- split_parallel(D, 1L)
    X <- split_parallel(X)

    ## Update with DTW in parallel
    D <- foreach(X = X,
                 distmat = D,
                 .combine = rbind,
                 .multicombine = TRUE,
                 .packages = "dtwclust",
                 .export = "enlist") %op% {
                     id_nn <- apply(distmat, 1L, which.min) # index of nearest neighbors
                     id_nn_prev <- id_nn + 1L # initialize all different
                     id_mat <- cbind(1L:nrow(distmat), id_nn) # to index the distance matrix

                     while (!is.null(y) && any(id_nn_prev != id_nn)) {
                         id_changed <- which(id_nn_prev != id_nn)
                         id_nn_prev <- id_nn

                         d_sub <- do.call(proxy::dist,
                                          enlist(x = X[id_changed],
                                                 y = Y[id_nn[id_changed]],
                                                 method = method,
                                                 dots = dots))

                         distmat[id_mat[id_changed, , drop = FALSE]] <- d_sub
                         id_nn <- apply(distmat, 1L, which.min)
                         id_mat[ , 2L] <- id_nn
                     }

                     distmat
                 }

    class(D) <- "crossdist"
    attr(D, "method") <- "DTW_LB"

    ## return
    D
}
