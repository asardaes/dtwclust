# ========================================================================================================
# Modified version of kcca to use lists of time series (to support different lengths)
# ========================================================================================================

kcca.list <- function (x, k, family = NULL, iter.max = 200L, trace = FALSE, ...)
{
     N <- length(x)

     if (N < k)
          stop("Number of clusters cannot be greater than number of observations in the data")

     centers <- x[sample(N,k)]
     cluster <- integer(N)
     k <- as.integer(k)

     for (iter in 1:iter.max) {
          clustold <- cluster
          distmat <- family@dist(x, centers, ...)
          cluster <- family@cluster(distmat = distmat)

          k <- length(centers)

          centers <- family@allcent(x, cluster = cluster, k = k, centers, ...)

          changes <- sum(cluster != clustold)

          if (trace) {
               td <- sum(distmat[cbind(1:N, cluster)])
               txt <- paste(changes, format(td), sep = " / ")
               cat("Iteration ", iter, ": ",
                   "Changes / Distsum = ",
                   formatC(txt, width = 12, format = "f"),
                   "\n", sep = "")
          }

          if (changes == 0) {
               if (trace) cat("\n")
               break
          }
     }

     cluster <- as.integer(family@cluster(distmat=distmat))
     cldist <- as.matrix(distmat[cbind(1:N, cluster)])
     size <- as.vector(table(cluster))
     clusinfo <- data.frame(size = size,
                            av_dist = as.vector(tapply(cldist[,1], cluster, sum))/size)

     list(cluster = cluster,
          centers = centers,
          clusinfo = clusinfo,
          cldist = cldist,
          iter = iter,
          converged = (iter < iter.max))
}
