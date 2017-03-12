# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

ddist2 <- function(distance, control) {
    symmetric <- if (is.null(control$symmetric)) FALSE else control$symmetric

    ## I need to re-register any custom distances in each parallel worker
    dist_entry <- proxy::pr_DB$get_entry(distance)

    ## Closures capture the values of the objects from the environment where they're created
    distfun <- function(x, centroids = NULL, ...) {
        if (!is.null(control$distmat)) {
            ## distmat pre-computed
            if (is.null(centroids)) {
                ## return whole distmat
                d <- control$distmat

            } else {
                ## distmat matrix already calculated, just subset it
                d <- control$distmat[ , attr(centroids, "id_cent"), drop = FALSE]
            }

        } else {
            ## distmat not available, calculate it
            ## Extra distance parameters in case of parallel computation
            ## They can be for the function or for proxy::dist
            dots <- list(...)

            ## Added defaults
            if (is.null(dots$window.size))
                dots$window.type <- "none"
            else if (is.null(dots$window.type))
                dots$window.type <- "slantedband"

            dots$error.check <- FALSE

            ## dtw uses L2 by default, but in dtwclust I want dtw to use L1 by default
            ## Important for multivariate series
            if (toupper(dist_entry$names[1L]) == "DTW" && is.null(dots$dist.method))
                dots$dist.method <- "L1"

            ## If the function doesn't have '...', remove invalid arguments from 'dots'
            valid_args <- names(dots)

            if (is.function(dist_entry$FUN)) {
                if (!has_dots(dist_entry$FUN))
                    valid_args <- union(names(formals(proxy::dist)), names(formals(dist_entry$FUN)))

            } else {
                valid_args <- names(formals(proxy::dist))
            }

            dots <- dots[intersect(names(dots), valid_args)]

            ## variables/functions from the parent environments that should be exported
            export <- c("distance", "dist_entry", "check_consistency", "enlist")

            if (is.null(centroids) && symmetric && !isTRUE(dots$pairwise)) {
                if (dist_entry$loop && foreach::getDoParWorkers() > 1L && isTRUE(foreach::getDoParName() != "doSEQ")) {
                    ## WHOLE SYMMETRIC DISTMAT WITH proxy LOOP IN PARALLEL
                    ## Only half of it is computed
                    ## I think proxy can do this if y = NULL, but not in parallel

                    ## strict pairwise as in proxy::dist doesn't make sense here, it's pairwise between pairs
                    dots$pairwise <- TRUE

                    pairs <- call_pairs(length(x), lower = FALSE)

                    pairs <- split_parallel(pairs, 1L)

                    d <- foreach(pairs = pairs,
                                 .combine = c,
                                 .multicombine = TRUE,
                                 .packages = control$packages,
                                 .export = export) %op% {

                                     if (!check_consistency(dist_entry$names[1L], "dist"))
                                         do.call(proxy::pr_DB$set_entry, dist_entry)

                                     ## 'dots' has all extra arguments that are valid
                                     dd <- do.call(proxy::dist,
                                                   enlist(x = x[pairs[ , 1L]],
                                                          y = x[pairs[ , 2L]],
                                                          method = distance,
                                                          dots = dots))

                                     dd
                                 }

                    rm("pairs")

                    D <- matrix(0, nrow = length(x), ncol = length(x))
                    D[upper.tri(D)] <- d
                    D <- t(D)
                    D[upper.tri(D)] <- d

                    d <- D
                    attr(d, "class") <- "crossdist"
                    attr(d, "dimnames") <- list(names(x), names(x))
                    rm("D")

                } else {
                    ## WHOLE SYMMETRIC DISTMAT WITH CUSTOM LOOP OR SEQUENTIAL proxy LOOP
                    ## maybe one of my distances, or one included in proxy by default, let it handle parallelization
                    d <- as.matrix(do.call(proxy::dist,
                                           enlist(x = x,
                                                  y = NULL,
                                                  method = distance,
                                                  dots = dots)))

                    class(d) <- "crossdist"
                }

            } else {
                ## WHOLE DISTMAT OR SUBDISTMAT OR NOT SYMMETRIC
                if (is.null(centroids))
                    centroids <- x

                dim_names <- list(names(x), names(centroids))

                x <- split_parallel(x)

                if (isTRUE(dots$pairwise)) {
                    centroids <- split_parallel(centroids)
                    validate_pairwise(x, centroids)
                    combine <- c

                } else {
                    centroids <- lapply(1L:foreach::getDoParWorkers(), function(dummy) centroids)
                    if (length(centroids) > length(x)) centroids <- centroids[1L:length(x)]
                    combine <- rbind
                }

                d <- foreach(x = x, centroids = centroids,
                             .combine = combine,
                             .multicombine = TRUE,
                             .packages = control$packages,
                             .export = export) %op% {
                                 if (!check_consistency(dist_entry$names[1L], "dist"))
                                     do.call(proxy::pr_DB$set_entry, dist_entry)

                                 ## 'dots' has all extra arguments that are valid
                                 dd <- do.call(proxy::dist,
                                               enlist(x = x,
                                                      y = centroids,
                                                      method = distance,
                                                      dots = dots))

                                 dd
                             }

                if (isTRUE(dots$pairwise)) {
                    attr(d, "class") <- "pairdist"

                } else {
                    attr(d, "class") <- "crossdist"
                    attr(d, "dimnames") <- dim_names
                }
            }
        }

        attr(d, "method") <- toupper(distance)
        attr(d, "call") <- NULL

        d
    }

    distfun
}
