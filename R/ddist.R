# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

ddist <- function(distance, control = new("dtwclustControl"), distmat = NULL) {
    needs_window <- c("dtw_lb", "lbk", "lbi")

    if (distance %in% needs_window)
        control@window.size <- check_consistency(control@window.size, "window")

    if (is.null(control@window.size))
        window.type <- "none"
    else
        window.type <- "slantedband"

    ## Closures capture the values of the objects from the environment where they're created
    distfun <- function(x, centroids = NULL, ...) {
        if (!is.null(distmat)) {
            if (is.null(centroids)) {
                d <- distmat

            } else {
                ## distmat matrix already calculated, just subset it
                id_XC <- attr(centroids, "id_cent")

                ## in case I messed up something when assigning id_cent
                # if (is.null(id_XC) || any(is.na(id_XC)))
                #     id_XC <- sapply(centroids, FUN = function(i.c) {
                #         i.row <- sapply(x, function(i.x) {
                #             if (length(i.x) == length(i.c))
                #                 ret <- all(i.x == i.c)
                #             else
                #                 ret <- FALSE
                #
                #             ret
                #         })
                #
                #         ## Take the first one in case a series is repeated in the dataset
                #         which(i.row)[1L]
                #     })

                d <- distmat[ , id_XC, drop = FALSE]
            }

        } else {
            ## distmat not available, calculate it
            ## Extra distance parameters in case of parallel computation
            ## They can be for the function or for proxy::dist
            dots <- list(...)

            dots$window.type <- window.type
            dots$window.size <- control@window.size
            dots$norm <- control@norm
            dots$error.check <- FALSE

            ## I need to re-register any custom distances in each parallel worker
            dist_entry <- proxy::pr_DB$get_entry(distance)

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
            export <- c("distance", "check_consistency", "enlist")

            if (is.null(centroids) && control@symmetric && !isTRUE(dots$pairwise)) {
                if (dist_entry$loop) {
                    ## WHOLE SYMMETRIC DISTMAT WITH proxy LOOP
                    ## Only half of it is computed
                    ## I think proxy can do this if y = NULL, but not in parallel

                    ## strict pairwise as in proxy::dist doesn't make sense here, but this case needs it
                    dots$pairwise <- TRUE

                    pairs <- call_pairs(length(x), lower = FALSE)

                    pairs <- split_parallel(pairs, 1L)

                    d <- foreach(pairs = pairs,
                                 .combine = c,
                                 .multicombine = TRUE,
                                 .packages = control@packages,
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
                    ## WHOLE SYMMETRIC DISTMAT WITH CUSTOM LOOP
                    ## most likely one of my distances, let it handle parallelization
                    d <- do.call(proxy::dist,
                                 enlist(x = x,
                                        y = NULL,
                                        method = distance,
                                        dots = dots))
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
                             .packages = control@packages,
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

                if (!is.null(dots$pairwise) && dots$pairwise) {
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
