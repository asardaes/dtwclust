# ========================================================================================================
# Modified version of flexclust::kcca to use lists of time series and/or support fuzzy clustering
# ========================================================================================================

kcca.list <- function (x, k, family, control, fuzzy = FALSE, ...) {
    N <- length(x)
    k <- as.integer(k)
    dots <- list(...)

    if (N < k)
        stop("Number of clusters cannot be greater than number of observations in the data")

    if (is.null(names(x)))
        names(x) <- paste0("series_", 1:N) # used by custom PAM centroids

    if (fuzzy) {
        cluster <- matrix(0, N, k)
        cluster[ , -1L] <- stats::runif(N *(k - 1)) / (k - 1)
        cluster[ , 1L] <- 1 - apply(cluster[ , -1L, drop = FALSE], 1L, sum)

        if (has_dots(family@allcent))
            centroids <- family@allcent(x = x, cl_id = cluster, k = k, cl_old = cluster, ...)
        else
            centroids <- do.call(family@allcent,
                                 enlist(x = x,
                                        cl_id = cluster,
                                        k = k,
                                        cl_old = cluster,
                                        dots = subset_dots(dots, family@allcent)))

    } else {
        id_cent <- sample(N, k)
        centroids <- x[id_cent]
        attr(centroids, "id_cent") <- id_cent
        cluster <- integer(N)
    }

    iter <- 1L
    objective_old <- Inf

    while (iter <= control@iter.max) {
        clustold <- cluster

        distmat <- family@dist(x, centroids, ...)

        cluster <- family@cluster(distmat = distmat, m = control@fuzziness)

        centroids <- family@allcent(x = x,
                                    cl_id = cluster,
                                    k = k,
                                    cent = centroids,
                                    cl_old = clustold,
                                    ...)

        if (fuzzy) {
            # fuzzy.R
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
        if (control@trace) cat("\n")
        warning("Clustering did not converge within the allowed iterations.")
        converged <- FALSE
        iter <- control@iter.max

    } else {
        converged <- TRUE
    }

    distmat <- family@dist(x, centroids, ...)
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
    size <- tabulate(cluster)

    ## if some clusters are empty, tapply() would not return enough rows
    clusinfo <- data.frame(size = size, av_dist = 0)
    clusinfo[clusinfo$size > 0L, "av_dist"] <- as.vector(tapply(cldist[ , 1L], cluster, mean))

    names(centroids) <- NULL
    attr(centroids, "id_cent") <- NULL
    centroids <- lapply(centroids, "attr<-", which = "id_cent", value = NULL)

    list(k = k,
         cluster = cluster,
         fcluster = fcluster,
         centroids = centroids,
         clusinfo = clusinfo,
         cldist = cldist,
         iter = iter,
         converged = converged)
}
