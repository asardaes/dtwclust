#' TADPole clustering
#'
#' Time-series Anytime Density Peaks Clustering as proposed by Begum et al. (2015).
#'
#' @export
#'
#' @param data A matrix or data frame where each row is a time series, or a list where each element
#'   is a time series. Multivariate series are not supported.
#' @param window.size Window size constraint for DTW (Sakoe-Chiba). See details.
#' @param k The number of desired clusters. Can be a vector with several values
#' @param dc The cutoff distance(s).
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
#' @template window
#'
#' @return A list with:
#'
#'   - `cl`: Cluster indices.
#'   - `centroids`: Indices of the centroids.
#'   - `distCalcPercentage`: Percentage of distance calculations that were actually performed.
#'
#' @template parallel
#'
#' @references
#'
#' Begum N, Ulanova L, Wang J and Keogh E (2015). ``Accelerating Dynamic Time Warping Clustering
#' with a Novel Admissible Pruning Strategy.'' In *Conference on Knowledge Discovery and Data
#' Mining*, series KDD '15. ISBN 978-1-4503-3664-2/15/08, \url{
#' http://dx.doi.org/10.1145/2783258.2783286}.
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
#' # Reinterpolate to same length
#' data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
#'
#' # Create parallel workers
#' cl <- makeCluster(detectCores())
#' invisible(clusterEvalQ(cl, library(dtwclust)))
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
#' # Compute VI Index
#' cat("VI index for TADPole:", cvi(kc.tadp$cl, CharTrajLabels, "VI"), "\n\n")
#' }
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

    ## Calculate matrices with bounds
    LBM <- proxy::dist(x, x,
                       method = lb,
                       window.size = window.size,
                       force.symmetry = TRUE,
                       norm = "L2",
                       error.check = error.check)

    ## NOTE: Euclidean is only valid as upper bound if 'symmetric1' step pattern is used
    UBM <- proxy::dist(x, x, method = "L2")

    `%this_op%` <- if (length(dc) >= foreach::getDoParWorkers()) `%op%` else foreach::`%do%`

    ## Return
    RET <- foreach(
        dc = dc,
        .combine = c,
        .multicombine = TRUE,
        .packages = c("dtwclust", "foreach"),
        .export = c("split_parallel", "%op%")) %this_op%
        {
            ## =====================================================================================
            ## Pruning during local density calculation
            ## =====================================================================================

            ## Flag definition
            # 0 - DTW calculated, and it lies below dc
            # 1 - calculate DTW
            # 2 - within dc, prune
            # 3 - not within dc, prune
            # 4 - identical series

            if (trace) cat("\tPruning during local density calculation\n")

            ## Initialize matrices
            D <- matrix(NA_real_, n, n)
            Flags <- matrix(-1L, n, n)

            ## Calculate values for upper triangular of the flag matrix (column-wise)
            utv <- unlist(lapply(2L:n, function(j) {
                sapply(1L:(j-1L), function(i) {
                    if (LBM[i,j] <= dc && UBM[i,j] > dc)
                        1L
                    else if (LBM[i,j] <= dc && UBM[i,j] < dc)
                        2L
                    else if (LBM[i,j] > dc)
                        3L
                    else
                        4L
                })
            }))

            ## NOTE: only upper triangular to avoid unnecessary DTW calculations
            Flags[upper.tri(Flags)] <- utv
            ## Calculate full DTW for those entries whose lower bounds lie around 'dc'
            ind1 <- which(Flags == 1L, arr.ind = TRUE)

            if (nrow(ind1) > 0L) {
                ind1 <- split_parallel(ind1, 1L)
                exclude <- setdiff(ls(), c("x", "window.size"))

                d1 <- foreach(ind1 = ind1,
                              .combine = c,
                              .multicombine = TRUE,
                              .packages = c("dtwclust", "foreach"),
                              .noexport = exclude) %op% {
                                  proxy::dist(x[ind1[ , 1L]], x[ind1[ , 2L]],
                                              method = "dtw_basic",
                                              window.size = window.size,
                                              step.pattern = symmetric1,
                                              norm = "L2",
                                              pairwise = TRUE)
                              }
                ## Fill distance matrix where necessary
                if (is.list(ind1)) ind1 <- do.call(rbind, ind1)
                D[ind1] <- d1
                Flags[ind1[d1 <= dc, ]] <- 0L # For 'Rho' calculation
            }

            ## Force symmetry
            D[lower.tri(D)] <- t(D)[lower.tri(D)]
            Flags[lower.tri(Flags)] <- t(Flags)[lower.tri(Flags)]

            ## Calculate local density vector
            Rho <- apply(Flags, 1L, function(i) sum(i == 0L | i == 2L))

            if (all(Rho == 0L))
                stop("No density peaks detected, choose a different value for cutoff distance 'dc'")

            if (max(Rho) == min(Rho))
                Rho <- rep(1L, length(Rho))
            else
                Rho <- (Rho - min(Rho)) / (max(Rho) - min(Rho))

            ## =====================================================================================
            ## Pruning during NN distance calculation from higher density list (phase 1)
            ## =====================================================================================

            if (trace) cat("\tPruning during nearest-neighbor distance calculation (phase 1)\n")

            RhoSorted <- sort(Rho, decreasing = TRUE, index.return = TRUE)

            deltaUB <- t(sapply(2L:n, function (i) {
                ## Index of higher density neighbors
                indHDN <- RhoSorted$ix[1L:(i-1L)]
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

            ## =====================================================================================
            ## Pruning during NN distance calculation from higher density list (phase 2)
            ## =====================================================================================

            if (trace) cat("\tPruning during nearest-neighbor distance calculation (phase 2)\n")

            # start at two
            i <- split_parallel(2L:n)

            DNN <- foreach(i = i,
                           .combine = rbind,
                           .multicombine = TRUE,
                           .packages = c("dtwclust", "foreach")) %op% {
                               t(sapply(i, function (i) {
                                   ## Index of higher density neighbors
                                   indHDN <- TADPorder$ix[1L:(i-1L)]
                                   ## Index of current object
                                   ii <- TADPorder$ix[i]

                                   ## If this is true, prune the calculation of true distance
                                   indPrune <- LBM[ii, indHDN] > deltaUB[ii]
                                   ## If the distance was already computed, don't do it again
                                   indPre <- Flags[ii, indHDN] == 0L | Flags[ii, indHDN] == 1L

                                   ## 'delta' will have the distances of HDN
                                   ## Initially filled with upper bound
                                   delta <- UBM[ii, indHDN]
                                   ## If some distances were already computed, put them here
                                   delta[indPre] <- D[ii, indHDN[indPre]]

                                   ## If the distance is not to be pruned nor previously calculated,
                                   ## compute it now
                                   indCompute <- !(indPrune | indPre)

                                   if (any(indCompute)) {
                                       d2 <- proxy::dist(x[ii], x[indHDN[indCompute]],
                                                         method = "dtw_basic",
                                                         window.size = window.size,
                                                         step.pattern = symmetric1,
                                                         norm = "L2")

                                       delta[indCompute] <- d2
                                   }

                                   c(delta = min(delta),
                                     NN = indHDN[which.min(delta)],
                                     distCalc = sum(indCompute))
                               }))
                           }

            ## To order according to default order
            indOrig <- sort(TADPorder$ix, index.return = TRUE)$ix
            delta <- c(max(DNN[ , 1L]), DNN[ , 1L])
            NN <- c(-1L, DNN[ , 2L])

            ## =====================================================================================
            ## Cluster assignment
            ## =====================================================================================

            if (trace) cat("\tPerforming cluster assignment\n\n")

            ## How many calculations were actually performed
            distCalc <-  sum(utv == 1L) + sum(DNN[ , 3L])
            distCalc <- (distCalc / (n * (n+1) / 2 - n)) * 100

            ## Normalize
            if (max(delta) == min(delta))
                zDelta <- rep(1L, length(delta))
            else
                zDelta <- (delta - min(delta)) / (max(delta) - min(delta))

            ret <- lapply(k, function(k) {
                ## Those with the most density are the cluster centroids (PAM)
                cent <- sort(sort(Rho * zDelta[indOrig],
                                  decreasing = TRUE,
                                  index.return = TRUE)$ix[1L:k])

                ## Assign a unique number to each cluster centroid
                cl <- rep(-1L, n)
                cl[cent] <- 1L:k

                ## Which elements don't have a label yet (this is ordered according to
                ## the density values of TADPorder, because the assignment is sequential)
                indCl <- TADPorder$ix

                ## Do the assignment (must be a loop)
                for (i in seq_along(indCl)) {
                    if (cl[indCl[i]] == -1L) cl[indCl[i]] <- cl[NN[i]]
                }

                if (any(cl == -1L))
                    warning(c("At least one series wasn't assigned to a cluster. ",
                              "This shouldn't happen, please contact maintainer."))

                if (trace)
                    cat("TADPole completed for k = ", k,
                        " & dc = ", dc,
                        ", pruning percentage = ",
                        formatC(100 - distCalc,
                                digits = 3L,
                                width = -1L,
                                format = "fg"),
                        "%\n\n",
                        sep = "")

                ## Return
                list(cl = cl,
                     centroids = cent,
                     distCalcPercentage = distCalc)
            })

            ## return from foreach
            ret
        }

    ## Return
    if (length(RET) == 1L) RET[[1L]] else RET
}
