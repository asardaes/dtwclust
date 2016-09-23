#' DTW Barycenter Averaging
#'
#' A global averaging method for time series under DTW (Petitjean, Ketterlin and Gancarski, 2011).
#'
#' This function tries to find the optimum average series between a group of time series in DTW space. Refer to
#' the cited article for specific details on the algorithm.
#'
#' If a given series reference is provided in \code{centroid}, the algorithm should always converge to the same
#' result provided the elements of \code{X} keep the same values, although their order may change.
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
#' Multivariate series should be provided as a list of matrices where time spans the rows and the variables
#' span the columns.
#' @param centroid Optionally, a time series to use as reference. Defaults to a random series of \code{X} if
#' \code{NULL}. For multivariate series, this should be a matrix with the same characteristics as the
#' matrices in \code{X}.
#' @param center Deprecated, please use \code{centroid} instead.
#' @param max.iter Maximum number of iterations allowed.
#' @param norm Norm for the local cost matrix of DTW. Either "L1" for Manhattan distance or "L2" for Euclidean
#' distance.
#' @param window.size Window constraint for the DTW calculations. \code{NULL} means no constraint. A slanted
#' band is used by default.
#' @param delta At iteration \code{i}, if \code{all(abs(centroid_{i}} \code{ - centroid_{i-1})} \code{ < delta)},
#' convergence is assumed.
#' @param error.check Should inconsistencies in the data be checked?
#' @param trace If \code{TRUE}, the current iteration is printed to screen.
#' @param ... Further arguments for \code{\link[dtw]{dtw}} or \code{\link{dtw_basic}}, e.g.
#' \code{step.pattern}.
#' @param dba.alignment Character indicating which function to use for calculating alignments, either
#' \code{\link[dtw]{dtw}} or \code{\link{dtw_basic}}. The latter should be faster.
#'
#' @return The average time series.
#'
#' @export
#'

DBA <- function(X, centroid = NULL, center = NULL, max.iter = 20L,
                norm = "L1", window.size = NULL, delta = 1e-3,
                error.check = TRUE, trace = FALSE, ..., dba.alignment = "dtw_basic") {
     if (!missing(center)) {
          warning("The 'center' argument has been deprecated, please use 'centroid' instead.")

          if (is.null(centroid)) centroid <- center
     }

     dba.alignment <- match.arg(dba.alignment, c("dtw", "dtw_basic"))

     X <- consistency_check(X, "tsmat")

     if (is.null(centroid))
          centroid <- X[[sample(length(X), 1L)]] # Random choice

     if (error.check) {
          consistency_check(X, "vltslist")
          consistency_check(centroid, "ts")
     }

     norm <- match.arg(norm, c("L1", "L2"))

     dots <- list(...)

     if (!is.null(window.size)) {
          w <- consistency_check(window.size, "window")

          if (is.null(dots$window.type))
               dots$window.type <- "slantedband"

     } else {
          w <- NULL
     }

     ## utils.R
     if (check_multivariate(X)) {
          ## multivariate
          mv <- reshape_multviariate(X, centroid) # utils.R

          new_c <- mapply(mv$series, mv$cent, SIMPLIFY = FALSE,
                          FUN = function(xx, cc) {
                               DBA(xx, cc,
                                   norm = norm,
                                   window.size = window.size,
                                   max.iter = max.iter,
                                   delta = delta,
                                   error.check = FALSE,
                                   trace = trace,
                                   dba.alignment = dba.alignment,
                                   ...)
                          })

          return(do.call(cbind, new_c))
     }

     ## pre-allocate local cost matrices
     if (dba.alignment == "dtw")
          LCM <- lapply(X, function(x) { matrix(0, length(x), length(centroid)) })
     else
          LCM <- lapply(X, function(x) { matrix(0, length(x) + 1L, length(centroid) + 1L) })

     ## maximum length of considered series
     L <- max(lengths(X))

     ## Register doSEQ if necessary
     check_parallel()

     Xs <- split_parallel(X)
     LCMs <- split_parallel(LCM)

     ## Iterations
     iter <- 1L
     centroid_old <- centroid

     if (trace) cat("\tDBA Iteration:")

     while(iter <= max.iter) {
          ## Return the coordinates of each series in X grouped by the coordinate they match to in the centroid time series
          ## Also return the number of coordinates used in each case (for averaging below)
          xg <- foreach(X = Xs, LCM = LCMs,
                        .combine = c,
                        .multicombine = TRUE,
                        .export = "enlist",
                        .packages = c("dtwclust", "stats")) %dopar% {
                             mapply(X, LCM, SIMPLIFY = FALSE, FUN = function(x, lcm) {
                                  if (dba.alignment == "dtw") {
                                       .Call("update_lcm", lcm, x, centroid,
                                             isTRUE(norm == "L2"), PACKAGE = "dtwclust")

                                       d <- do.call(dtw::dtw, enlist(x = lcm, window.size = w, dots = dots))

                                  } else {
                                       d <- do.call(dtw_basic, enlist(x = x, y = centroid,
                                                                      window.size = w, norm = norm,
                                                                      backtrack = TRUE, gcm = lcm,
                                                                      dots = dots))
                                  }

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

          if (all(abs(centroid - centroid_old) < delta)) {
               if (trace)
                    cat("", iter ,"- Converged!\n")

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
