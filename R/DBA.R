#' DTW Barycenter Averaging
#'
#' A global averaging method for time series under DTW (Petitjean, Ketterlin and Gancarski, 2011).
#'
#' This function tries to find the optimum average series between a group of time series in DTW space. Refer to
#' the cited article for specific details on the algorithm.
#'
#' If a given series reference is provided in \code{center}, the algorithm should always converge to the same
#' result provided the rows of \code{X} keep the same values, although their order may change.
#'
#' @references
#'
#' Petitjean F, Ketterlin A and Gancarski P (2011). ``A global averaging method for dynamic time warping, with applications to
#' clustering.'' \emph{Pattern Recognition}, \strong{44}(3), pp. 678 - 693. ISSN 0031-3203, \url{
#' http://dx.doi.org/10.1016/j.patcog.2010.09.013}, \url{
#' http://www.sciencedirect.com/science/article/pii/S003132031000453X}.
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
#' all(dtw.avg == dtw.avg2)
#'
#' @param X A data matrix where each row is a time series. Optionally, a list where each element is a time series.
#' @param center Optionally, a time series to use as reference. It must be a numeric vector. Defaults to a
#' random series of \code{X} if \code{NULL}.
#' @param max.iter Maximum number of iterations allowed.
#' @param error.check Should inconsistencies in the data be checked?
#' @param trace If \code{TRUE}, the current iteration is printed to screen.
#'
#' @return The average time series.
#'
#' @export
#' @importFrom stats aggregate
#' @importFrom dtw dtw

DBA <- function(X, center = NULL, max.iter = 50, error.check = TRUE, trace = FALSE) {

     X <- consistency_check(X, "tsmat")

     n <- length(X)

     if (is.null(center))
          center <- X[[sample(n, 1)]] # Random choice

     if (error.check) {
          consistency_check(X, "vltslist")
          consistency_check(center, "ts")
     }

     iter <- 1
     C.old <- center-1

     while(iter<=max.iter) {

          ## Return the coordinates of each series in X grouped by the coordinate they match to in the center time series
          ## Also return the number of coordinates used in each case (for averaging below)
          xg <- lapply(X, function(x) {
               d <- dtw::dtw(x, center)

               x.sub <- stats::aggregate(x[d$index1], by=list(ind = d$index2), sum)

               n.sub <- stats::aggregate(x[d$index1], by=list(ind = d$index2), length)

               cbind(sum = x.sub$x, n = n.sub$x)
          })

          ## Put everything in one big data frame
          xg <- reshape2::melt(xg)

          ## Aggregate according to index of center time series (Var1) and also the variable type (Var2)
          xg <- stats::aggregate(xg$value, by = list(xg$Var1, xg$Var2), sum)

          ## Average
          center <- xg$x[xg$Group.2 == "sum"] / xg$x[xg$Group.2 == "n"]

          if (all(center == C.old)) {
               iter <- iter+1

               if (trace)
                    cat("DBA: Iteration", iter ,"- Converged!\n\n")

               break

          } else {
               iter <- iter+1
               C.old <- center

               if (trace)
                    cat("DBA: Iteration", iter, "\n")
          }
     }

     if (iter > max.iter)
          warning("DBA algorithm did not converge within the allowed iterations.")

     as.numeric(center)
}
