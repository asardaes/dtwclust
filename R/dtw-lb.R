#' DTW distance matrix guided by Lemire's improved lower bound
#'
#' Calculation of a distance matrix with the Dynamic Time Warping (DTW) distance guided by Lemire's
#' improved lower bound (LB_Improved).
#'
#' This function first calculates an initial estimate of a distance matrix between two sets of time series
#' using LB_Improved. Afterwards, it uses the estimate to calculate the corresponding true DTW distance
#' between \emph{only} the nearest neighbors of each series in \code{x} found in \code{y}, and it continues
#' iteratively until no changes in the nearest neighbors occur.
#'
#' If only \code{x} is provided, the distance matrix is calculated between all its time series.
#'
#' This could be useful in case one is interested in only the nearest neighbor of one or more series
#' within a dataset.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10(2) + 1 = 21} observations falling within
#' the window.
#'
#' @section Parallel Computing:
#'
#' Please note that running tasks in parallel does \strong{not} guarantee faster computations.
#' The overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#' The user can register a parallel backend, e.g. with the \code{doParallel} package, in order to attempt to
#' speed up the calculations (see the examples).
#'
#' @note
#'
#' This function uses a lower bound that is only defined for time series of equal length.
#'
#' @seealso
#'
#' \code{\link{lb_keogh}}, \code{\link{lb_improved}}
#'
#' @references
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .''
#' \emph{Pattern Recognition}, \strong{42}(9), pp. 2169 - 2180. ISSN 0031-3203,
#' \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030},
#' \url{http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
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
#' NN1 <- apply(d, 1, which.min)
#'
#' # Calculate the DTW distances between all elements (about three times slower)
#' system.time(d2 <- proxy::dist(data[1:5], data[6:50], method = "DTW",
#'                               window.type = "slantedband", window.size = 20))
#'
#' # Nearest neighbors
#' NN2 <- apply(d2, 1, which.min)
#'
#' # Same results?
#' all(NN1 == NN2)
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
#' @author Alexis Sarda-Espinosa
#'
#' @param x,y A matrix where rows are time series, or a list of time series.
#' @param window.size Window size to use with the LB and DTW calculation. See details.
#' @param norm Pointwise distance. Either \code{"L1"} for Manhattan distance or \code{"L2"} for Euclidean.
#' @param error.check Should inconsistencies in the data be checked?
#' @param pairwise Calculate pairwise distances?
#' @param dtw.func Which function to use for core DTW the calculations, either "dtw" or "dtw_basic". See
#' \code{\link[dtw]{dtw}} and \code{\link{dtw_basic}}.
#' @param ... Further arguments for \code{dtw.func} or \code{\link{lb_improved}}.
#'
#' @return The distance matrix with class \code{crossdist}.
#'
#' @export

dtw_lb <- function(x, y = NULL, window.size = NULL, norm = "L1",
                   error.check = TRUE, pairwise = FALSE,
                   dtw.func = "dtw_basic", ...) {
     norm <- match.arg(norm, c("L1", "L2"))
     dtw.func <- match.arg(dtw.func, c("dtw", "dtw_basic"))

     if (dtw.func == "dtw")
          method <- ifelse(norm == "L1", "DTW", "DTW2")
     else
          method <- toupper(dtw.func)

     X <- consistency_check(x, "tsmat")

     check_parallel()

     dots <- list(...)
     dots$dist.method <- norm
     dots$norm <- norm
     dots$window.size <- window.size
     dots$pairwise <- TRUE

     if (pairwise) {
          if (is.null(y))
               Y <- x
          else
               Y <- consistency_check(y, "tsmat")

          if (length(X) != length(Y))
               stop("Pairwise distances require the same amount of series in 'x' and 'y'")

          if (is.null(window.size))
               dots$window.type <- "none"
          else
               dots$window.type <- "slantedband"

          X <- split_parallel(X)
          Y <- split_parallel(Y)

          D <- foreach(X = X, Y = Y,
                       .combine = c,
                       .multicombine = TRUE,
                       .packages = "dtwclust",
                       .export = "enlist") %dopar% {
                            do.call(proxy::dist,
                                    enlist(x = X, y = Y,
                                           method = method,
                                           dots = dots))
                       }

          return(D)
     }

     window.size <- consistency_check(window.size, "window")
     dots$window.size <- window.size
     dots$window.type <- "slantedband"

     ## NOTE: I tried starting with LBK estimate, refining with LBI and then DTW but, overall,
     ## it was usually slower, almost the whole matrix had to be recomputed for LBI

     if (!is.null(y))
          Y <- consistency_check(y, "tsmat")
     else
          Y <- X

     ## Initial estimate
     D <- proxy::dist(X, Y, method = "LBI",
                      window.size = window.size,
                      norm = norm,
                      error.check = error.check,
                      ...)

     ## Update with DTW
     new.indNN <- apply(D, 1L, which.min) # index of nearest neighbors
     indNN <- new.indNN + 1L # initialize all different

     while (any(new.indNN != indNN)) {
          indNew <- which(new.indNN != indNN)
          indNN <- new.indNN

          indNew <- split_parallel(indNew)

          exclude <- setdiff(ls(), c("X", "Y", "method", "dots", "indNew", "indNN"))

          dSub <- foreach(indNew = indNew,
                          .combine = c,
                          .multicombine = TRUE,
                          .packages = "dtwclust",
                          .noexport = exclude,
                          .export = "enlist") %dopar% {
                               do.call(proxy::dist,
                                       enlist(x = X[indNew],
                                              y = Y[indNN[indNew]],
                                              method = method,
                                              dots = dots))
                          }

          D[cbind(1L:length(X), indNN)[unlist(indNew), , drop = FALSE]] <- dSub

          new.indNN <- apply(D, 1L, which.min)
     }

     class(D) <- "crossdist"
     attr(D, "method") <- "DTW_LB"

     D
}
