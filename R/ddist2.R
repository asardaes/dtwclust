# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

ddist2 <- function(distance, control) {
    symmetric <- isTRUE(control$symmetric)
    ## I need to re-register any custom distances in each parallel worker
    dist_entry <- proxy::pr_DB$get_entry(distance)

    ## Closures capture the values of the objects from the environment where they're created
    distfun <- function(x, centroids = NULL, ...) {
        x <- any2list(x)
        if (!is.null(centroids)) centroids <- any2list(centroids)

        if (!is.null(control$distmat)) {
            if (inherits(control$distmat, "Distmat")) {
                ## internal class, sparse or full
                i <- 1L:length(x)
                j <- if (is.null(centroids)) i else control$distmat$id_cent
                d <- control$distmat[i, j, drop = FALSE]

            } else {
                ## normal distmat pre-computed
                if (is.null(centroids))
                    d <- control$distmat ## return whole distmat
                else
                    d <- control$distmat[, control$distmat$id_cent, drop = FALSE] ## subset
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
            export <- c("distance", "dist_entry", "check_consistency", "enlist", "subset_dots")

            if (is.null(centroids) && symmetric && !isTRUE(dots$pairwise)) {
                if (dist_entry$loop && foreach::getDoParWorkers() > 1L) {
                    ## WHOLE SYMMETRIC DISTMAT WITH proxy LOOP IN PARALLEL
                    ## Only half of it is computed
                    ## I think proxy can do this if y = NULL, but not in parallel
                    len <- length(x)
                    loop_endpoints <- symmetric_loop_endpoints(len)
                    d <- bigmemory::big.matrix(len, len, "double", 0)
                    d_desc <- bigmemory::describe(d)

                    foreach(
                        loop_endpoints = loop_endpoints,
                        .combine = c,
                        .multicombine = TRUE,
                        .noexport = c("d"),
                        .packages = c(control$packages, "bigmemory"),
                        .export = export
                    ) %op% {
                        dd <- bigmemory::attach.big.matrix(d_desc)
                        ## while should be faster here, no big data is modified
                        i <- loop_endpoints$start$i
                        j <- loop_endpoints$start$j
                        while (j <= loop_endpoints$end$j) {
                            i_max <- if (j == loop_endpoints$end$j) loop_endpoints$end$i else len
                            while (i <= i_max) {
                                dist_val <- do.call(dist_entry$FUN,
                                                    enlist(x[[i]],
                                                           x[[j]],
                                                           dots = subset_dots(dots,
                                                                              dist_entry$FUN)))
                                dd[i,j] <- dist_val
                                dd[j,i] <- dist_val
                                i <- i + 1L
                            }
                            j <- j + 1L
                            i <- j + 1L
                        }
                        rm("dd")
                        gc()
                        NULL
                    }

                    d <- d[,]
                    attr(d, "class") <- "crossdist"
                    attr(d, "dimnames") <- list(names(x), names(x))
                    gc()

                } else {
                    ## WHOLE SYMMETRIC DISTMAT WITH CUSTOM LOOP OR SEQUENTIAL proxy LOOP
                    ## maybe one of my distances, or one included in proxy by default, let it handle parallelization
                    d <- base::as.matrix(do.call(proxy::dist,
                                                 enlist(x = x,
                                                        y = NULL,
                                                        method = distance,
                                                        dots = dots)))

                    class(d) <- "crossdist"
                }

            } else {
                ## WHOLE DISTMAT OR SUBDISTMAT OR NOT SYMMETRIC
                if (is.null(centroids)) centroids <- x
                dim_names <- list(names(x), names(centroids))
                x <- split_parallel(x)

                if (isTRUE(dots$pairwise)) {
                    centroids <- split_parallel(centroids)
                    validate_pairwise(x, centroids)
                    combine <- c

                } else {
                    centroids <- lapply(1L:foreach::getDoParWorkers(), function(dummy) { centroids })
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

        ## return
        d
    }

    distfun
}
