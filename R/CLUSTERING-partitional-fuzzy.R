# ==================================================================================================
# Fuzzy objective function
# ==================================================================================================

fuzzy_objective <- function(u, distmat, m) { sum(u^m * distmat^2) }

# ==================================================================================================
# Partitional/fuzzy clustering
# ==================================================================================================

#' @importFrom stats runif
#'
pfclust <- function (x, k, family, control, fuzzy = FALSE, cent, trace = FALSE, args) {
    N <- length(x)
    k <- as.integer(k)
    if (is.null(control$version)) control$version <- 2L # nocov

    if (fuzzy && cent == "fcm") {
        cluster <- matrix(0, N, k)
        cluster[, -1L] <- stats::runif(N * (k - 1L)) / (k - 1L)
        cluster[, 1L] <- 1 - apply(cluster[, -1L, drop = FALSE], 1L, sum)
        centroids <- quoted_call(
            family@allcent,
            x = x,
            cl_id = cluster,
            k = k,
            cl_old = cluster,
            dots = subset_dots(args$cent, family@allcent)
        )
    }
    else {
        id_cent <- sample(N, k)
        centroids <- x[id_cent]
        if (inherits(control$distmat, "Distmat")) control$distmat$id_cent <- id_cent
        cluster <- integer(N)
    }

    iter <- 1L
    clustold <- cluster
    objective_old <- Inf
    dmi <- cbind(1L:N, integer(N))
    while (iter <= control$iter.max) {
        distmat <- quoted_call(family@dist, x = x, centroids = centroids, dots = args$dist)
        cluster <- family@cluster(distmat = distmat, m = control$fuzziness)
        if (control$version < 2L) { # nocov start
            centroids <- quoted_call(
                family@allcent,
                x = x,
                cl_id = cluster,
                k = k,
                cent = centroids,
                cl_old = clustold,
                dots = subset_dots(args$cent, family@allcent)
            )
        } # nocov end

        # NOTE: a custom fuzzy centroid function that doesn't want to rely on the usual fuzzy
        # objective would not be able to customize this convergence criterion...
        if (fuzzy && cent != "fcmdd") {
            objective <- fuzzy_objective(cluster, distmat = distmat, m = control$fuzziness)
            if (trace) {
                cat("Iteration ", iter, ": ",
                    "Objective = ",
                    formatC(objective, width = 6, format = "f"),
                    "\n", sep = "")
            }
            if (abs(objective - objective_old) < control$delta) {
                if (trace) cat("\n")
                break
            }
            objective_old <- objective
        }
        else {
            dmi[, 2L] <- if (cent != "fcmdd") cluster else apply(cluster, 1L, which.max)
            changes <- sum(dmi[, 2L] != clustold)
            if (trace) {
                td <- sum(distmat[dmi])
                txt <- paste(changes, format(td), sep = " / ")
                cat("Iteration ", iter, ": ",
                    "Changes / Distsum = ",
                    formatC(txt, width = 12, format = "f"),
                    "\n", sep = "")
            }
            if (changes == 0L) {
                if (trace) cat("\n")
                break
            }
        }
        iter <- iter + 1L
        if (control$version > 1L)
            centroids <- quoted_call(
                family@allcent,
                x = x,
                cl_id = cluster,
                k = k,
                cent = centroids,
                cl_old = clustold,
                dots = subset_dots(args$cent, family@allcent)
            )

        clustold <- dmi[, 2L]
    }

    if (iter > control$iter.max) {
        if (trace) cat("\n") # nocov
        converged <- FALSE
        iter <- control$iter.max
    }
    else {
        converged <- TRUE
    }

    if (control$version < 2L) { # nocov start
        distmat <- quoted_call(family@dist, x = x, centroids = centroids, dots = args$dist)
        cluster <- family@cluster(distmat = distmat, m = control$fuzziness)
    } # nocov end

    if (fuzzy) {
        fcluster <- cluster
        cluster <- max.col(-distmat, "first")
        rownames(fcluster) <- names(x)
        colnames(fcluster) <- paste0("cluster_", 1L:k)
    }
    else {
        fcluster <- matrix(NA_real_)
    }

    cldist <- base::as.matrix(distmat[cbind(1L:N, cluster)])
    clusinfo <- compute_clusinfo(k, cluster, cldist) # UTILS-utils.R

    names(centroids) <- NULL
    if (inherits(control$distmat, "Distmat")) {
        attr(centroids, "series_id") <- base::unname(control$distmat$id_cent)
    }

    # return
    list(k = k,
         cluster = cluster,
         fcluster = fcluster,
         centroids = centroids,
         clusinfo = clusinfo,
         cldist = cldist,
         iter = iter,
         converged = converged)
}
