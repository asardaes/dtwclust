#' TADPole clustering
#'
#' Time-series Anytime Density Peaks Clustering as proposed by Begum et al., 2015.
#'
#' This function can be called either directly or through \code{\link{dtwclust}}.
#'
#' TADPole clustering adopts a relatively new clustering framework and adapts it to time series clustering
#' with DTW. See the cited article for the details of the algorithm.
#'
#' Because of the way the algorithm works, it can be considered a kind of Partitioning Around Medoids (PAM).
#' This means that the cluster centers are always elements of the data.
#'
#' The algorithm first uses the DTW's upper and lower bounds to find series with many close neighbors (in
#' DTW space). Anything below the cutoff distance (\code{dc}) is considered a neighbor. Aided with this
#' information, the algorithm then tries to prune as many DTW calculations as possible in order to accelerate
#' the clustering procedure. The series that lie in dense areas (i.e. that have lots of neighbors) are taken
#' as cluster centers.
#'
#' The algorithm relies on the DTW bounds, which are only defined for time series of equal lengths.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10*2 + 1 = 21} observations falling within
#' the window.
#'
#' @section Parallel Computing:
#'
#' Please note that running tasks in parallel does \strong{not} guarantee faster computations.
#' The overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#' The user can register a parallel backend with the \code{doParallel} package in order to attempt to
#' speed up the calculations (see the examples).
#'
#' @examples
#'
#' \dontrun{
#' #### Running TADPole with parallel support
#' require(doParallel)
#'
#' # Load data
#' data(uciCT)
#'
#' # Reinterpolate to same length and coerce as matrix
#' data <- t(sapply(CharTraj, reinterpolate, newLength = 180))
#'
#' # Create parallel workers
#' cl <- makeCluster(detectCores())
#' registerDoParallel(cl)
#'
#' # Cluster
#' kc.tadp <- TADPole(data, k = 20, window.size = 20, dc = 1.5)
#'
#' # Stop parallel workers
#' stopCluster(cl)
#'
#' # Return to sequential computations
#' registerDoSEQ()
#'
#' # Compute Rand Index
#' cat("Rand index for TADPole:", randIndex(kc.tadp$cl, CharTrajLabels), "\n\n")
#' }
#'
#' @references
#'
#' Begum N, Ulanova L, Wang J and Keogh E (2015). ``Accelerating Dynamic Time Warping Clustering with a Novel Admissible Pruning
#' Strategy.'' In \emph{Conference on Knowledge Discovery and Data Mining}, series KDD '15. ISBN 978-1-4503-3664-2/15/08, \url{
#' http://dx.doi.org/10.1145/2783258.2783286}.
#'
#' @param data The data matrix where each row is a time series. Optionally a list with each time series.
#' @param window.size Window size constraint for DTW. See details.
#' @param k The number of desired clusters.
#' @param dc The cutoff distance.
#' @param error.check Should the data be checked for inconsistencies?
#'
#' @return A list with: \itemize{
#'   \item \code{cl}: Cluster indices.
#'   \item \code{centers}: Indices of the centers.
#'   \item \code{distCalcPercentage}: Percentage of distance calculations that were actually performed.
#' }
#'
#' @import foreach
#' @importFrom parallel splitIndices
#'
#' @export

TADPole <- function(data, window.size = NULL, k = 2, dc, error.check = TRUE) {

     if (missing(window.size))
          stop("Please provide a positive window size")
     if (missing(dc))
          stop("Please provide the 'dc' parameter")

     x <- consistency_check(data, "tsmat")

     n <- length(x)

     if (n < 2)
          stop("data should have more than one time series")
     if (k > n)
          stop("Number of clusters should be less than the number of time series")

     ## Calculate matrices with bounds
     LBM <- proxy::dist(x, x, method = "LBK",
                        window.size = window.size, force.symmetry = TRUE,
                        norm = "L2",
                        error.check = error.check)

     UBM <- proxy::dist(x, x, method = "L2")

     ## Euclidean is only valid as upper bound if 'symmetric1' step pattern is used
     step.pattern <- get("symmetric1") # so that CHECK doesn't complain

     ## Attempt parallel computations?
     do_par <- foreach::getDoParRegistered()

     if (do_par)
          workers <- foreach::getDoParWorkers()
     else
          workers <- 1L

     ## ============================================================================================================================
     ## Pruning during local density calculation
     ## ============================================================================================================================

     ## Flag definition
     # 0 - DTW calculated, and it lies below dc
     # 1 - calculate DTW
     # 2 - within dc, prune
     # 3 - not within dc, prune
     # 4 - identical series

     ## Initialize matrices
     D <- matrix(NA, n, n)
     Flags <- matrix(-1L, n, n)

     ## Calculate values for the upper triangular part of the flag matrix (column-wise)
     utv <- unlist(lapply(2:n, function(j) {

          sapply(1:(j-1), function(i) {

               if (LBM[i,j]<=dc && UBM[i,j]>dc)
                    f <- 1L
               else if (LBM[i,j]<=dc && UBM[i,j]<dc)
                    f <- 2L
               else if (LBM[i,j] > dc)
                    f <- 3L
               else
                    f <- 4L

               f
          })
     }))

     ## NOTE: only the upper triangular is filled here to avoid unnecessary calculations of DTW
     Flags[upper.tri(Flags)] <- utv

     ## Calculate full DTW for those entries whose lower bounds lie around 'dc'
     ind1 <- which(Flags == 1, arr.ind = TRUE)

     if (nrow(ind1) > 0) {

          ## Attempt parallel calculations?
          if (do_par) {
               tasks <- parallel::splitIndices(nrow(ind1), workers)

               ind1 <- lapply(tasks, function(id) {
                    ind1[id,]
               })

               exclude <- setdiff(ls(), c("x", "step.pattern", "window.size"))

               d1 <- foreach(ind1 = ind1,
                             .combine = c,
                             .multicombine = TRUE,
                             .packages = "dtwclust",
                             .noexport = exclude) %dopar% {
                                  proxy::dist(x[ind1[, 1]], x[ind1[, 2]], method = "DTW2",
                                              step.pattern = step.pattern,
                                              window.type = "slantedband", window.size = window.size,
                                              pairwise = TRUE)
                             }
          } else {
               d1 <- proxy::dist(x[ind1[, 1]], x[ind1[, 2]], method = "DTW2",
                                 step.pattern = step.pattern,
                                 window.type = "slantedband", window.size = window.size,
                                 pairwise = TRUE)
          }

          ## Fill distance matrix where necessary
          D[which(Flags==1)] <- d1
          Flags[which(Flags==1)[d1 <= dc]] <- 0 # For 'Rho' calculation
     }

     ## Force symmetry
     D[lower.tri(D)] <- t(D)[lower.tri(D)]
     Flags[lower.tri(Flags)] <- t(Flags)[lower.tri(Flags)]

     ## Calculate local density vector
     Rho <- apply(Flags, 1, function(i) sum(i==0 | i==2))

     if (all(Rho==0))
          stop("No density peaks detected, choose a different value for cutoff distance 'dc'")

     if (max(Rho) == min(Rho))
          Rho <- rep(1L, length(Rho))
     else
          Rho <- (Rho - min(Rho)) / (max(Rho) - min(Rho))

     ## ============================================================================================================================
     ## Pruning during NN distance calculation from higher density list (phase 1)
     ## ============================================================================================================================

     RhoSorted <- sort(Rho, decreasing = TRUE, index.return = TRUE)

     deltaUB <- t(sapply(2:n, function (i) {
          ## Index of higher density neighbors
          indHDN <- RhoSorted$ix[1:(i-1)]
          ## Index of current object
          ii <- RhoSorted$ix[i]

          ## Use existing distance if available, otherwise use the upper bound
          ub.i <- D[ii, indHDN]
          indDNA <- is.na(ub.i)
          ub.i[indDNA] <- UBM[ii, indHDN[indDNA]]

          min(ub.i)
     }))

     deltaUB <- c(max(deltaUB), deltaUB)

     ## Reorder the values
     indOrig <- sort(RhoSorted$ix, index.return = TRUE)$ix
     deltaUB <- deltaUB[indOrig]

     ## New order
     TADPorder <- sort(Rho * deltaUB, decreasing = TRUE, index.return = TRUE)

     ## ============================================================================================================================
     ## Pruning during NN distance calculation from higher density list (phase 2)
     ## ============================================================================================================================

     ## Attempt parallel calculations?
     if (do_par) {
          tasks <- parallel::splitIndices(n-1L, workers)

          i <- lapply(tasks, function(id) {
               id + 1 # start at two
          })

          DNN <- foreach(i = i,
                         .combine = rbind,
                         .multicombine = TRUE,
                         .packages = "dtwclust") %dopar% {

                              t(sapply(i, function (i) {
                                   ## Index of higher density neighbors
                                   indHDN <- TADPorder$ix[1:(i-1)]
                                   ## Index of current object
                                   ii <- TADPorder$ix[i]

                                   ## If this is true, prune the calculation of the true distance
                                   indPrune <- LBM[ii, indHDN] > deltaUB[ii]
                                   ## If the distance was already computed before, don't do it again
                                   indPre <- Flags[ii, indHDN] == 0 | Flags[ii, indHDN] == 1

                                   ## 'delta' will have the distances to neighbors with higher densities. Initially filled with upper bound
                                   delta <- UBM[ii, indHDN]
                                   ## If some distances were already computed, put them here
                                   delta[indPre] <- D[ii, indHDN[indPre]]

                                   ## If the distance is not to be pruned nor previously calculated, compute it now
                                   indCompute <- !(indPrune | indPre)

                                   if (sum(indCompute) > 0) {
                                        d2 <- proxy::dist(x[ii], x[indHDN[indCompute]], method = "DTW2",
                                                          step.pattern = step.pattern,
                                                          window.type = "slantedband", window.size = window.size)

                                        delta[indCompute] <- d2
                                   }

                                   c(delta = min(delta), NN = indHDN[which.min(delta)], distCalc = sum(indCompute))
                              }))
                         }

     } else {
          DNN <- t(sapply(2:n, function (i) {
               ## Index of higher density neighbors
               indHDN <- TADPorder$ix[1:(i-1)]
               ## Index of current object
               ii <- TADPorder$ix[i]

               ## If this is true, prune the calculation of the true distance
               indPrune <- LBM[ii, indHDN] > deltaUB[ii]
               ## If the distance was already computed before, don't do it again
               indPre <- Flags[ii, indHDN] == 0 | Flags[ii, indHDN] == 1

               ## 'delta' will have the distances to neighbors with higher densities. Initially filled with upper bound
               delta <- UBM[ii, indHDN]
               ## If some distances were already computed, put them here
               delta[indPre] <- D[ii, indHDN[indPre]]

               ## If the distance is not to be pruned nor previously calculated, compute it now
               indCompute <- !(indPrune | indPre)

               if (sum(indCompute) > 0) {
                    ## TODO: parallel support
                    d2 <- proxy::dist(x[ii], x[indHDN[indCompute]], method = "DTW2",
                                      step.pattern = step.pattern,
                                      window.type = "slantedband", window.size = window.size)

                    delta[indCompute] <- d2
               }

               c(delta = min(delta), NN = indHDN[which.min(delta)], distCalc = sum(indCompute))
          }))
     }

     ## To order according to default order
     indOrig <- sort(TADPorder$ix, index.return = TRUE)$ix

     delta <- c(max(DNN[, 1]), DNN[, 1])

     NN <- c(-1, DNN[, 2])

     ## ============================================================================================================================
     ## Cluster assignment
     ## ============================================================================================================================

     ## Normalize
     if(max(delta) == min(delta))
          zDelta <- rep(1L, length(delta))
     else
          zDelta <- (delta - min(delta)) / (max(delta) - min(delta))

     ## Those with the most density are the cluster centers (PAM)
     C <- sort(Rho * zDelta[indOrig], decreasing = TRUE, index.return = TRUE)$ix[1:k]
     C <- sort(C)

     ## Assign a unique number to each cluster center
     cl <- rep(-1L, n)
     cl[C] <- 1L:k

     ## Which elements don't have a label yet (this is ordered according to the density values of TADPorder,
     ## because the assignment is sequential)
     indCl <- TADPorder$ix

     ## Do the assignment (must be a loop)
     for (i in seq_along(indCl)) {
          if (cl[indCl[i]] == -1)
               cl[indCl[i]] <- cl[NN[i]]
     }

     ## How many calculations were actually performed
     distCalc <-  sum(utv == 1) + sum(DNN[,3])

     if(any(cl == -1))
          warning(c("At least one series wasn't assigned to a cluster. ",
                    "This shouldn't happen, please contact maintainer."))

     ## Return
     list(cl = cl,
          centers = C,
          distCalcPercentage = (distCalc / (n * (n+1) / 2 - n)) * 100)
}
