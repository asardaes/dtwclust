# ========================================================================================================
# Modified version of kcca to use lists of time series (to support different lengths)
# ========================================================================================================

kcca.list <- function (x, k, family = NULL, iter.max = 30L, trace = FALSE, ...)
{
     N <- length(x)

     if (N < k)
          stop("Number of clusters cannot be greater than number of observations in the data")

     id_cen <- sample(N,k)
     centers <- x[id_cen]
     #attr(centers, "id") <- id_cen
     cluster <- integer(N)
     k <- as.integer(k)
     iter <- 1L

     while (iter <= iter.max) {
          clustold <- cluster
          distmat <- family@dist(x, centers, ...)
          cluster <- family@cluster(distmat = distmat)

          centers <- family@allcent(x, cluster, k, centers, clustold, ...)

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

          iter <- iter + 1L
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
          converged = (iter <= iter.max))
}
