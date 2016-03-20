# ========================================================================================================
# Modified version of kcca to use lists of time series (to support different lengths)
# ========================================================================================================

kcca.list <- function (x, k, family, control, fuzzy = FALSE, ...)
{
     N <- length(x)
     k <- as.integer(k)

     if (N < k)
          stop("Number of clusters cannot be greater than number of observations in the data")

     if (is.null(names(x)))
          names(x) <- paste0("series_", 1:N) # used by custom PAM centers

     if (fuzzy) {
          cluster <- matrix(0, N, k)
          cluster[ , -1L] <- stats::runif(N *(k - 1L)) / (k - 1)
          cluster[ , 1L] <- 1 - apply(cluster[ , -1L], 1L, sum)
          centers <- family@allcent(x, cluster, k, ...)

     } else {
          id_cent <- sample(N, k)
          centers <- x[id_cent]
          attr(centers, "id_cent") <- id_cent
          cluster <- integer(N)
     }

     iter <- 1L
     objective_old <- Inf

     while (iter <= control@iter.max) {
          clustold <- cluster

          distmat <- family@dist(x, centers, ...)

          cluster <- family@cluster(distmat = distmat, m = control@fuzziness)

          centers <- family@allcent(x, cluster, k, centers, clustold, ...)

          if (fuzzy) {
               # utils.R
               objective <- fuzzy_objective(cluster, distmat = distmat, m = control@fuzziness)

               if (control@trace) {
                    cat("Iteration ", iter, ": ",
                        "Objective = ",
                        formatC(objective, width = 6, format = "f"),
                        "\n", sep = "")
               }

               if (abs(objective - objective_old) < control@delta) {
                    if (control@trace) cat("\n")
                    break
               }

               objective_old <- objective

          } else {
               changes <- sum(cluster != clustold)

               if (control@trace) {
                    td <- sum(distmat[cbind(1L:N, cluster)])
                    txt <- paste(changes, format(td), sep = " / ")
                    cat("Iteration ", iter, ": ",
                        "Changes / Distsum = ",
                        formatC(txt, width = 12, format = "f"),
                        "\n", sep = "")
               }

               if (changes == 0L) {
                    if (control@trace) cat("\n")
                    break
               }
          }

          iter <- iter + 1L
     }

     if (iter > control@iter.max) {
          warning("Clustering did not converge within the allowed iterations.")
          converged <- FALSE
          iter <- control@iter.max

     } else {
          converged <- TRUE
     }

     distmat <- family@dist(x, centers, ...)
     cluster <- family@cluster(distmat = distmat, m = control@fuzziness)

     if (fuzzy) {
          fcluster <- cluster
          rownames(fcluster) <- names(x)
          colnames(fcluster) <- paste0("cluster_", 1:k)

          cluster <- max.col(-distmat, "first")

     } else {
          fcluster <- matrix(NA_real_)
     }

     cldist <- as.matrix(distmat[cbind(1L:N, cluster)])
     size <- as.vector(table(cluster))
     clusinfo <- data.frame(size = size,
                            av_dist = as.vector(tapply(cldist[ , 1L], cluster, sum))/size)

     names(centers) <- NULL
     attr(centers, "id_cent") <- NULL
     centers <- lapply(centers, "attr<-", which = "id_cent", value = NULL)

     list(k = k,
          cluster = cluster,
          fcluster = fcluster,
          centers = centers,
          clusinfo = clusinfo,
          cldist = cldist,
          iter = iter,
          converged = converged)
}
