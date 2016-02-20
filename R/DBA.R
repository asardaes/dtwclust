#' DTW Barycenter Averaging
#'
#' A global averaging method for time series under DTW (Petitjean, Ketterlin and Gancarski, 2011).
#'
#' This function tries to find the optimum average series between a group of time series in DTW space. Refer to
#' the cited article for specific details on the algorithm.
#'
#' If a given series reference is provided in \code{center}, the algorithm should always converge to the same
#' result provided the elements of \code{X} keep the same values, although their order may change.
#'
#' @section Parallel Computing:
#'
#' Please note that running tasks in parallel does \strong{not} guarantee faster computations.
#' The overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#' The user can register a parallel backend, e.g. with the \code{doParallel} package, in order to attempt to
#' speed up the calculations (see the examples).
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
#' plot(dtw.avg, type="l")
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
#' @param X A data matrix where each row is a time series, or a list where each element is a time series.
#' @param center Optionally, a time series to use as reference. It must be a numeric vector. Defaults to a
#' random series of \code{X} if \code{NULL}.
#' @param max.iter Maximum number of iterations allowed.
#' @param norm Norm for the local cost matrix of DTW. Either "L1" for Manhattan distance or "L2" for Euclidean
#' distance.
#' @param window.size Window constraint for the DTW calculations. \code{NULL} means no constraint.
#' @param delta At iteration \code{i}, if \code{all(abs(center_{i}} \code{ - center_{i-1})} \code{ < delta)},
#' convergence is assumed.
#' @param error.check Should inconsistencies in the data be checked?
#' @param trace If \code{TRUE}, the current iteration is printed to screen.
#' @param ... Further arguments for \code{\link[dtw]{dtw}}, e.g. \code{step.pattern}. Do not provide
#' \code{window.type} here, just set \code{window.size} to the desired value.
#'
#' @return The average time series.
#'
#' @export
#'

DBA <- function(X, center = NULL, max.iter = 20L,
                norm = "L1", window.size = NULL, delta = 1e-3,
                error.check = TRUE, trace = FALSE, ...) {

     X <- consistency_check(X, "tsmat")

     if (!is.null(window.size)) {
          w <- consistency_check(window.size, "window")
          window.type = "slantedband"

     } else {
          w <- NULL
          window.type = "none"
     }

     norm <- match.arg(norm, c("L1", "L2"))

     ## for C helper
     square <- norm == "L2"

     n <- length(X)

     if (is.null(center))
          center <- X[[sample(n, 1L)]] # Random choice

     if (error.check) {
          consistency_check(X, "vltslist")
          consistency_check(center, "ts")
     }

     ## pre-allocate local cost matrices
     LCM <- lapply(X, function(x) { matrix(0, length(x), length(center)) })

     ## maximum length of considered series
     L <- max(lengths(X))

     ## Register doSEQ if necessary
     check_parallel()

     Xs <- split_parallel(X)
     LCMs <- split_parallel(LCM)

     dots <- list(...)

     ## Iterations
     iter <- 1L
     center_old <- center

     while(iter <= max.iter) {
          ## Return the coordinates of each series in X grouped by the coordinate they match to in the center time series
          ## Also return the number of coordinates used in each case (for averaging below)
          xg <- foreach(X = Xs, LCM = LCMs,
                        .combine = c,
                        .multicombine = TRUE,
                        .packages = c("dtwclust", "stats")) %dopar% {
                             mapply(X, LCM, SIMPLIFY = FALSE, FUN = function(x, lcm) {
                                  .Call("update_lcm", lcm, x, center, square, PACKAGE = "dtwclust")

                                  d <- do.call(dtw::dtw, c(list(x = lcm,
                                                                window.type = window.type,
                                                                window.size = w),
                                                           dots))

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

          ## Aggregate according to index of center time series (Var1) and also the variable type (Var2)
          xg <- stats::aggregate(xg$value, by = list(xg$Var1, xg$Var2), sum)

          ## Average
          center <- xg$x[xg$Group.2 == "sum"] / xg$x[xg$Group.2 == "n"]

          if (all(abs(center - center_old) < delta)) {
               if (trace)
                    cat("DBA: Iteration", iter ,"- Converged!\n\n")

               break

          } else {
               center_old <- center

               if (trace)
                    cat("DBA: Iteration", iter, "\n")

               iter <- iter + 1L
          }
     }

     if (iter > max.iter)
          warning("DBA algorithm did not converge within the allowed iterations.")

     as.numeric(center)
}
