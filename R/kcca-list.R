# ========================================================================================================
# Modified version of kcca to use lists of time series (to support different lengths)
# ========================================================================================================

kcca.list <- function (x, k, family = kccaFamily("kmeans"), control = NULL)
{
     MYCALL <- match.call()
     control <- as(control, "flexclustControl")
     x <- family@preproc(x)
     N <- length(x)

     if (k < 2)
          stop("Number of clusters must be at least 2")
     if (N < k)
          stop("Number of clusters cannot be greater than number of observations in the data")

     if (control@classify == "auto") {
          control@classify = "hard"
     }

     centers <- x[sample(N,k)]
     cluster <- integer(N)
     k <- as.integer(k)

     sannprob <- control@simann[1]

     if (control@classify %in% c("hard", "simann")) {
          for (iter in 1:control@iter.max) {
               clustold <- cluster
               distmat <- family@dist(x, centers)
               cluster <- family@cluster(x, distmat = distmat)

               if (control@classify == "simann") {
                    #cluster <- flexclust:::perturbClusters(cluster, sannprob)
                    #sannprob <- sannprob * control@simann[2]
                    stop("Unimplemented here")
               }

               k <- ifelse(is.matrix(centers), nrow(centers), length(centers))

               centers <- family@allcent(x, cluster = cluster, k = k)

               changes <- sum(cluster != clustold)

               if (control@verbose && (iter%%control@verbose == 0)) {
                    td <- sum(distmat[cbind(1:N, cluster)])
                    txt <- paste(changes, format(td), sep = " / ")
                    cat(formatC(iter, width = 6),
                        "Changes / Distsum :",
                        formatC(txt, width = 12, format = "f"),
                        "\n")
               }

               if (changes == 0)
                    break
          }
     }
     else if (control@classify == "weighted") {
          stop("Unsupported for series of different lengths")
     }
     else
          stop("Unknown classification method")

     cluster <- as.integer(family@cluster(distmat=distmat))
     cldist <- as.matrix(distmat[cbind(1:N, cluster)])
     size <- as.vector(table(cluster))
     clusinfo <- data.frame(size = size,
                            av_dist = as.vector(tapply(cldist[,1], cluster, sum))/size)

     z <- new("kccasimple",
              k = as.integer(k),
              cluster = cluster,
              centers = centers,
              family = family,
              clusinfo = clusinfo,
              cldist = cldist,
              iter = iter,
              converged = (iter < control@iter.max),
              call = MYCALL,
              control = control)

     z
}
