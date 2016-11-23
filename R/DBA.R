#' DTW Barycenter Averaging
#'
#' A global averaging method for time series under DTW (Petitjean, Ketterlin and Gancarski, 2011).
#'
#' @export
#'
#' @param X A data matrix where each row is a time series, or a list where each element is a time
#'   series. Multivariate series should be provided as a list of matrices where time spans the rows
#'   and the variables span the columns.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of
#'   \code{X} if \code{NULL}. For multivariate series, this should be a matrix with the same
#'   characteristics as the matrices in \code{X}.
#' @param ... Currently ignored.
#' @param window.size Window constraint for the DTW calculations. \code{NULL} means no constraint. A
#'   slanted band is used by default.
#' @param norm Norm for the local cost matrix of DTW. Either "L1" for Manhattan distance or "L2" for
#'   Euclidean distance.
#' @param max.iter Maximum number of iterations allowed.
#' @param delta At iteration \code{i}, if \code{all(abs(centroid_{i}} \code{ - centroid_{i-1})}
#'   \code{ < delta)}, convergence is assumed.
#' @param error.check Should inconsistencies in the data be checked?
#' @param trace If \code{TRUE}, the current iteration is printed to screen.
#'
#' @details
#'
#' This function tries to find the optimum average series between a group of time series in DTW
#' space. Refer to the cited article for specific details on the algorithm.
#'
#' If a given series reference is provided in \code{centroid}, the algorithm should always converge
#' to the same result provided the elements of \code{X} keep the same values, although their order
#' may change.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in
#' \code{window.size} that represents the distance between the point considered and one of the edges
#' of the window. Therefore, if, for example, \code{window.size = 10}, the warping for an
#' observation \eqn{x_i} considers the points between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting
#' in \code{10(2) + 1 = 21} observations falling within the window.
#'
#' @return The average time series.
#'
#' @section Parallel Computing:
#'
#'   Please note that running tasks in parallel does \strong{not} guarantee faster computations. The
#'   overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#'   The user can register a parallel backend, e.g. with the \code{doParallel} package, in order to
#'   attempt to speed up the calculations (see the examples).
#'
#' @references
#'
#' Petitjean F, Ketterlin A and Gancarski P (2011). ``A global averaging method for dynamic time
#' warping, with applications to clustering.'' \emph{Pattern Recognition}, \strong{44}(3), pp. 678 -
#' 693. ISSN 0031-3203, \url{http://dx.doi.org/10.1016/j.patcog.2010.09.013},
#' \url{http://www.sciencedirect.com/science/article/pii/S003132031000453X}.
#'
#' @examples
#'
#' # Sample data
#' data(uciCT)
#'
#' # Obtain an average for the first 5 time series
#' dtw.avg <- DBA(CharTraj[1:5], CharTraj[[1]], trace = TRUE)
#'
#' # Plot
#' matplot(do.call(cbind, CharTraj[1:5]), type = "l")
#' points(dtw.avg)
#'
#' # Change the provided order
#' dtw.avg2 <- DBA(CharTraj[5:1], CharTraj[[1]], trace = TRUE)
#'
#' # Same result?
#' all(dtw.avg == dtw.avg2)
#'
#' \dontrun{
#' #### Running DBA with parallel support
#' # For such a small dataset, this is probably slower in parallel
#' require(doParallel)
#'
#' # Create parallel workers
#' cl <- makeCluster(detectCores())
#' invisible(clusterEvalQ(cl, library(dtwclust)))
#' registerDoParallel(cl)
#'
#' # DTW Average
#' cen <- DBA(CharTraj[1:5], CharTraj[[1]], trace = TRUE)
#'
#' # Stop parallel workers
#' stopCluster(cl)
#'
#' # Return to sequential computations
#' registerDoSEQ()
#' }
#'
DBA <- function(X, centroid = NULL, ...,
                window.size = NULL, norm = "L1",
                max.iter = 20L, delta = 1e-3,
                error.check = TRUE, trace = FALSE)
{
    X <- any2list(X)

    if (is.null(centroid))
        centroid <- X[[sample(length(X), 1L)]] # Random choice

    if (error.check) {
        check_consistency(X, "vltslist")
        check_consistency(centroid, "ts")
    }

    ## utils.R
    if (is_multivariate(X)) {
        ## multivariate
        mv <- reshape_multviariate(X, centroid) # utils.R

        new_c <- mapply(mv$series, mv$cent, SIMPLIFY = FALSE,
                        FUN = function(xx, cc) {
                            DBA(xx, cc, ...,
                                norm = norm,
                                window.size = window.size,
                                max.iter = max.iter,
                                delta = delta,
                                error.check = FALSE,
                                trace = trace)
                        })

        return(do.call(cbind, new_c))
    }

    norm <- match.arg(norm, c("L1", "L2"))

    if (!is.null(window.size))
        window.size <- check_consistency(window.size, "window")

    dots <- list(...)

    ## maximum length of considered series
    L <- max(lengths(X))

    Xs <- split_parallel(X)

    ## pre-allocate local cost matrices
    GCM <- NULL # for CHECK
    GCMs <- lapply(Xs, function(dummy) { list(matrix(0, L + 1L, length(centroid) + 1L)) })

    ## Iterations
    iter <- 1L
    centroid_old <- centroid

    if (trace) cat("\tDBA Iteration:")

    while(iter <= max.iter) {
        ## Return the coordinates of each series in X grouped by the coordinate they match to in the centroid time series
        ## Also return the number of coordinates used in each case (for averaging below)
        xg <- foreach(X = Xs, GCM = GCMs,
                      .combine = c,
                      .multicombine = TRUE,
                      .export = "enlist",
                      .packages = c("dtwclust", "stats")) %op% {
                          mapply(X, GCM, SIMPLIFY = FALSE, FUN = function(x, gcm) {
                              d <- do.call(dtw_basic,
                                           enlist(x = x, y = centroid,
                                                  window.size = window.size, norm = norm,
                                                  backtrack = TRUE, gcm = gcm,
                                                  dots = dots))

                              x.sub <- stats::aggregate(x[d$index1],
                                                        by = list(ind = d$index2),
                                                        sum)

                              n.sub <- stats::aggregate(x[d$index1],
                                                        by = list(ind = d$index2),
                                                        length)

                              cbind(sum = x.sub$x, n = n.sub$x)
                          })
                      }

        ## Put everything in one big data frame
        xg <- reshape2::melt(xg)

        ## Aggregate according to index of centroid time series (Var1) and also the variable type (Var2)
        xg <- stats::aggregate(xg$value, by = list(xg$Var1, xg$Var2), sum)

        ## Average
        centroid <- xg$x[xg$Group.2 == "sum"] / xg$x[xg$Group.2 == "n"]

        if (isTRUE(all.equal(abs(centroid), abs(centroid_old), tolerance = delta))) {
            if (trace) cat("", iter ,"- Converged!\n")

            break

        } else {
            centroid_old <- centroid

            if (trace) {
                cat(" ", iter, ",", sep = "")
                if (iter %% 10 == 0) cat("\n\t\t")
            }

            iter <- iter + 1L
        }
    }

    if (iter > max.iter && trace)
        cat(" Did not 'converge'\n")

    as.numeric(centroid)
}
