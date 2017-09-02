#' DTW distance matrix guided by Lemire's improved lower bound
#'
#' Calculation of a distance matrix with the Dynamic Time Warping (DTW) distance guided by Lemire's
#' improved lower bound (LB_Improved).
#'
#' @export
#'
#' @param x,y A matrix or data frame where rows are time series, or a list of time series.
#' @param window.size Window size to use with the LB and DTW calculation. See details.
#' @param norm Either `"L1"` for Manhattan distance or `"L2"` for Euclidean.
#' @template error-check
#' @param pairwise Calculate pairwise distances?
#' @param dtw.func Which function to use for core DTW the calculations, either "dtw" or "dtw_basic".
#'   See [dtw::dtw()] and [dtw_basic()].
#' @param nn.margin Either 1 to search for nearest neighbors row-wise, or 2 to search column-wise.
#'   Only implemented for `dtw.func` = "dtw_basic".
#' @param ... Further arguments for `dtw.func` or [lb_improved()].
#'
#' @details
#'
#' This function first calculates an initial estimate of a distance matrix between two sets of time
#' series using [lb_improved()] (the [proxy::dist()] version). Afterwards, it uses the estimate to
#' calculate the corresponding true DTW distance between *only* the nearest neighbors of each series
#' in `x` found in `y`, and it continues iteratively until no changes in the nearest neighbors
#' occur.
#'
#' If only `x` is provided, the distance matrix is calculated between all its time series,
#' effectively returning a matrix filled with the LB_Improved values.
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
#' with [proxy::dist()]. Also note that [lb_improved()] calculates warping envelopes for the series
#' in `y`, so be careful with the provided order and `nn.margin` (see examples).
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
#' system.time(d <- dtw_lb(data[1:5], data[6:50], window.size = 20L))
#'
#' # Nearest neighbors
#' NN1 <- apply(d, 1L, which.min)
#'
#' # Calculate the DTW distances between all elements (slower)
#' system.time(d2 <- proxy::dist(data[1:5], data[6:50], method = "DTW",
#'                               window.type = "sakoechiba", window.size = 20L))
#'
#' # Nearest neighbors
#' NN2 <- apply(d2, 1L, which.min)
#'
#' # Calculate the DTW distances between all elements using dtw_basic
#' # (actually faster, see notes)
#' system.time(d3 <- proxy::dist(data[1:5], data[6:50], method = "DTW_BASIC",
#'                               window.size = 20L))
#'
#' # Nearest neighbors
#' NN3 <- apply(d3, 1L, which.min)
#'
#' # Change order and margin for nearest neighbor search
#' # (usually fastest, see notes)
#' system.time(d4 <- dtw_lb(data[6:50], data[1:5],
#'                          window.size = 20L, nn.margin = 2L))
#'
#' # Nearest neighbors *column-wise*
#' NN4 <- apply(d4, 2L, which.min)
#'
#' # Same results?
#' identical(NN1, NN2)
#' identical(NN1, NN3)
#' identical(NN1, NN4)
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
                   dtw.func = "dtw_basic", nn.margin = 1L, ...)
{
    norm <- match.arg(norm, c("L1", "L2"))
    dtw.func <- match.arg(dtw.func, c("dtw", "dtw_basic"))
    if (nn.margin != 1L) nn.margin <- 2L

    if (dtw.func == "dtw")
        method <- if (norm == "L1") "DTW" else "DTW2"
    else
        method <- toupper(dtw.func)

    X <- tslist(x)
    Y <- if (is.null(y)) X else tslist(y)

    if (is_multivariate(X) || is_multivariate(Y))
        stop("dtw_lb does not support multivariate series.")

    dots <- list(...)
    dots$dist.method <- norm
    dots$norm <- norm
    dots$window.size <- window.size
    dots$pairwise <- TRUE

    if (pairwise) {
        if (error.check) {
            check_consistency(X, "tslist")
            check_consistency(Y, "tslist")
        }

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
                                        dots = dots),
                                 TRUE)
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

    ## y = NULL means diagonal is zero and NNs are themselves
    if (is.null(y)) {
        class(D) <- "crossdist"
        attr(D, "method") <- "DTW_LB"
        return(D)
    }

    D <- split_parallel(D, 1L)
    X <- split_parallel(X)

    ## Update with DTW in parallel
    D <- foreach(X = X,
                 distmat = D,
                 .combine = rbind,
                 .multicombine = TRUE,
                 .packages = "dtwclust",
                 .export = c("enlist", "call_dtwlb")) %op% {
                     if (method == "DTW_BASIC") {
                         ## modifies distmat in place
                         dots$margin <- nn.margin
                         do.call(call_dtwlb,
                                 enlist(x = X, y = Y, distmat = distmat, dots = dots),
                                 TRUE)

                     } else {
                         if (nn.margin != 1L)
                             warning("Column-wise nearest neighbors are not implemented for dtw::dtw")
                         id_nn <- apply(distmat, 1L, which.min) # index of nearest neighbors
                         id_nn_prev <- id_nn + 1L # initialize all different
                         id_mat <- cbind(1L:nrow(distmat), id_nn) # to index the distance matrix
                         while (any(id_nn_prev != id_nn)) {
                             id_changed <- which(id_nn_prev != id_nn)
                             id_nn_prev <- id_nn

                             d_sub <- do.call(proxy::dist,
                                              enlist(x = X[id_changed],
                                                     y = Y[id_nn[id_changed]],
                                                     method = method,
                                                     dots = dots),
                                              TRUE)

                             distmat[id_mat[id_changed, , drop = FALSE]] <- d_sub
                             id_nn <- apply(distmat, 1L, which.min)
                             id_mat[ , 2L] <- id_nn
                         }
                     }

                     ## return from foreach()
                     distmat
                 }

    class(D) <- "crossdist"
    attr(D, "method") <- "DTW_LB"

    ## return
    D
}

## helper function
call_dtwlb <- function(x, y, distmat, ..., window.size, norm, margin,
                       step.pattern = NULL, gcm = NULL)
{
    if (is.null(step.pattern) || identical(step.pattern, symmetric2))
        step.pattern <- 2
    else if (identical(step.pattern, symmetric1))
        step.pattern <- 1
    else
        stop("step.pattern must be either symmetric1 or symmetric2 (without quotes)")

    L <- max(lengths(y)) + 1L

    if (is.null(gcm))
        gcm <- matrix(0, 2L, L)
    else if (!is.matrix(gcm) || nrow(gcm) < 2L || ncol(gcm) < L)
        stop("dtw_lb: Dimension inconsistency in 'gcm'")
    else if (storage.mode(gcm) != "double")
        stop("dtw_lb: If provided, 'gcm' must have 'double' storage mode.")

    dots <- list(window.size = window.size,
                 norm = switch(norm, "L1" = 1, "L2" = 2),
                 step.pattern = step.pattern,
                 backtrack = FALSE,
                 gcm = gcm)

    ## return
    .Call(C_dtw_lb, x, y, distmat, margin, dots, PACKAGE = "dtwclust")
}
