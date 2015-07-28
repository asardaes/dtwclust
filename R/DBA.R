#' DTW Barycenter Averaging
#'
#' See Petitjean 2011
#'
#' @param X A data matrix where each row is a time series. Optionally, a list where each element is a time series.
#' @param center Optionally, a time series to use as reference. Defaults to a random series of \code{X} if \code{NULL}.
#' @param max.iter Maximum number of iterations allowed.
#' @param error.check Should inconsistencies in the data be checked?
#' @param trace If \code{TRUE}, the current iteration is printed to screen.
#'
#' @return The average time series.
#'
#' @export

DBA <- function(X, center = NULL, max.iter = 25, error.check = TRUE, trace = FALSE) {

     ## For looping convenience
     if (is.matrix(X))
          X <- lapply(seq_len(nrow(X)), function(i) X[i,])

     n <- length(X)

     if (is.null(center))
          center <- X[[sample(n, 1)]] # Random choice

     if (error.check) {
          consistency_check(X, "tslist")
          consistency_check(center, "ts")
     }

     iter <- 0
     C.old <- center-1

     while(iter<=max.iter) {

          ## Return the coordinates of each series in X grouped by the coordinate they match to in the center time series
          ## Also return the number of coordinates used in each case (for averaging below)
          xg <- lapply(X, function(x) {
               d <- dtw(x, center)

               x.sub <- aggregate(x[d$index1], by=list(ind = d$index2), sum)

               n.sub <- aggregate(x[d$index1], by=list(ind = d$index2), length)

               cbind(ind = x.sub$ind, sum = x.sub$x, n = n.sub$x)
          })

          ## Put everything in one big data frame
          xg <- melt(xg) # from reshape2

          ## Aggregate according to index of center time series (Var1) and also the variable type (Var2)
          xg <- aggregate(xg$value, by = list(xg$Var1, xg$Var2), sum)

          ## Average
          center <- xg$x[xg$Group.2 == "sum"] / xg$x[xg$Group.2 == "n"]

          if (all(center == C.old)) {
               iter <- iter+1

               if (trace)
                    cat("Iteration", iter ,"- Converged!\n\n")

               break

          } else {
               iter <- iter+1
               C.old <- center

               if (trace)
                    cat("Iteration", iter, "\n")
          }
     }

     center
}
