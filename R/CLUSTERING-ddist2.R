# ==================================================================================================
# Helpers
# ==================================================================================================

# Use existing distmat if available
use_distmat <- function(distmat, x, centroids) {
    if (!inherits(distmat, "Distmat"))
        stop("Invalid distance matrix in control.") # nocov
    # internal class, sparse or full
    i <- 1L:length(x)
    j <- if (is.null(centroids)) i else distmat$id_cent
    distmat[i, j, drop = FALSE]
}

# Extra distance parameters in case of parallel computation
# They can be for the function or for proxy::dist
#' @importFrom proxy dist
#'
get_dots <- function(dist_entry, x, centroids, ...) {
    dots <- list(...)

    # Added defaults
    if (is.null(dots$window.size)) {
        dots$window.type <- "none"
    }
    else if (is.null(dots$window.type)) {
        dots$window.type <- "slantedband"
    }

    dots$error.check <- FALSE

    # dtw uses L2 by default, but in dtwclust I want dtw to use L1 by default
    # Important for multivariate series
    if (tolower(dist_entry$names[1L]) == "dtw" && is.null(dots$dist.method) && is_multivariate(c(x, centroids))) {
        dots$dist.method <- "L1" # nocov
    }

    # If the function doesn't have '...', remove invalid arguments from 'dots'
    valid_args <- names(dots)
    if (is.function(dist_entry$FUN)) {
        if (!has_dots(dist_entry$FUN)) {
            valid_args <- union(names(formals(proxy::dist)), names(formals(dist_entry$FUN)))
        }
    }
    else {
        valid_args <- names(formals(proxy::dist))
    }

    dots[intersect(names(dots), valid_args)]
}

# Function to split indices for the symmetric, parallel, proxy case
#' @importFrom parallel splitIndices
#'
split_parallel_symmetric <- function(n, num_workers, adjust = 0L) {
    if (num_workers <= 2L || n <= 4L) {
        mid_point <- as.integer(n / 2)
        # indices for upper part of the lower triangular
        ul_trimat <- 1L:mid_point + adjust
        # indices for lower part of the lower triangular
        ll_trimat <- (mid_point + 1L):n + adjust
        # put triangular parts together for load balance
        trimat <- list(ul = ul_trimat, ll = ll_trimat)
        attr(trimat, "trimat") <- TRUE
        trimat <- list(trimat)
        mid_point <- mid_point + adjust
        attr(ul_trimat, "rows") <- ll_trimat
        mat <- list(ul_trimat)
        ids <- c(trimat, mat)
    }
    else {
        mid_point <- as.integer(n / 2)
        # recursion
        rec1 <- split_parallel_symmetric(mid_point, as.integer(num_workers / 4), adjust)
        rec2 <- split_parallel_symmetric(n - mid_point, as.integer(num_workers / 4), mid_point + adjust)
        endpoints <- parallel::splitIndices(mid_point, max(length(rec1) + length(rec2), num_workers))
        endpoints <- endpoints[lengths(endpoints) > 0L]
        mat <- lapply(endpoints, function(ids) {
            ids <- ids + adjust
            attr(ids, "rows") <- (mid_point + 1L):n + adjust
            ids
        })
        ids <- c(rec1, rec2, mat)
    }
    chunk_sizes <- unlist(lapply(ids, function(x) {
        if (is.null(attr(x, "trimat"))) length(x) else median(lengths(x))
    }))
    # return
    ids[sort(chunk_sizes, index.return = TRUE)$ix]
}

# calculate only half the distance matrix in parallel
#' @importFrom proxy dist
#'
parallel_symmetric <- function(d_desc, ids, x, distance, dots) {
    attach.big.matrix <- get("attach.big.matrix", asNamespace("bigmemory"), mode = "function")
    dd <- attach.big.matrix(d_desc)
    if (isTRUE(attr(ids, "trimat"))) {
        # assign upper part of lower triangular
        ul <- ids$ul
        if (length(ul) > 1L) {
            dd[ul,ul] <- base::as.matrix(quoted_call(
                proxy::dist,
                x = x[ul],
                y = NULL,
                method = distance,
                dots = dots
            ))
        }

        # assign lower part of lower triangular
        ll <- ids$ll
        if (length(ll) > 1L) {
            dd[ll,ll] <- base::as.matrix(quoted_call(
                proxy::dist,
                x = x[ll],
                y = NULL,
                method = distance,
                dots = dots
            ))
        }
    }
    else {
        # assign matrix chunks
        rows <- attr(ids, "rows")
        mat_chunk <- base::as.matrix(quoted_call(
            proxy::dist,
            x = x[rows],
            y = x[ids],
            method = distance,
            dots = dots
        ))

        dd[rows,ids] <- mat_chunk
        dd[ids,rows] <- t(mat_chunk)
    }
}

# ==================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ==================================================================================================

#' @importFrom proxy dist
#' @importFrom proxy pr_DB
#'
ddist2 <- function(distance, control) {
    # I need to re-register any custom distances in each parallel worker
    dist_entry <- proxy::pr_DB$get_entry(distance)
    symmetric <- isTRUE(control$symmetric)
    warned <- FALSE

    # variables/functions from the parent environments that should be exported
    export <- c("check_consistency", "do_call", "quoted_call", "parallel_symmetric", "distance", "dist_entry")

    ret <- function(result, ...) {
        ret <- structure(result, method = toupper(distance), ...)
        if (!is.null(attr(ret, "call"))) {
            attr(ret, "call") <- NULL
        }
        ret
    }

    # Closures capture the values of the objects from the environment where they're created
    distfun <- function(x, centroids = NULL, ...) {
        x <- tslist(x)
        if (!is.null(centroids)) centroids <- tslist(centroids)

        if (length(x) == 1L && is.null(centroids)) {
            return(ret(base::matrix(0, 1L, 1L),
                       class = "crossdist",
                       dimnames = list(names(x), names(x))))
        }

        if (!is.null(control$distmat)) {
            return(ret(use_distmat(control$distmat, x, centroids)))
        }

        dots <- get_dots(dist_entry, x, centroids, ...)

        if (!dist_entry$loop) {
            # CUSTOM LOOP, LET THEM HANDLE OPTIMIZATIONS
            dm <- quoted_call(
                proxy::dist, x = x, y = centroids, method = distance, dots = dots
            )

            if (isTRUE(dots$pairwise)) {
                dim(dm) <- NULL
                return(ret(dm, class = "pairdist"))
            }
            else if (inherits(dm, "dist")) {
                return(ret(dm))
            }
            else {
                return(ret(base::as.matrix(dm), class = "crossdist"))
            }
        }

        if (is.null(centroids) && symmetric && !isTRUE(dots$pairwise)) {
            multiple_workers <- foreach::getDoParWorkers() > 1L

            if (multiple_workers && isNamespaceLoaded("bigmemory")) {
                # WHOLE SYMMETRIC DISTMAT IN PARALLEL
                # Only half of it is computed
                # proxy can do this if y = NULL, but not in parallel
                len <- length(x)

                # undo bigmemory's seed change, backwards reproducibility
                seed <- get0(".Random.seed", .GlobalEnv, mode = "integer")
                big.matrix <- get("big.matrix", asNamespace("bigmemory"), mode = "function")
                bigmemory_describe <- get("describe", asNamespace("bigmemory"), mode = "function")
                d <- big.matrix(len, len, "double", 0)
                d_desc <- bigmemory_describe(d)
                assign(".Random.seed", seed, .GlobalEnv)

                ids <- integer() # "initialize", so CHECK doesn't complain about globals
                foreach(
                    ids = split_parallel_symmetric(len, foreach::getDoParWorkers()),
                    .combine = c,
                    .multicombine = TRUE,
                    .noexport = c("d"),
                    .packages = c(control$packages, "bigmemory"),
                    .export = export
                ) %op% {
                    if (!check_consistency(dist_entry$names[1L], "dist")) {
                        do_call(proxy::pr_DB$set_entry, dist_entry) # nocov
                    }

                    parallel_symmetric(d_desc, ids, x, distance, dots)
                    NULL
                }

                # coerce to normal matrix
                return(ret(d[,], class = "crossdist", dimnames = list(names(x), names(x))))
            }
            else if (multiple_workers && !warned && isTRUE(getOption("dtwclust_suggest_bigmemory", TRUE))) {
                warned <<- TRUE
                warning("Package 'bigmemory' is not available, cannot parallelize computation with '",
                        distance,
                        "'. Use options(dtwclust_suggest_bigmemory = FALSE) to avoid this warning.")
            }
            else if (!multiple_workers) {
                # WHOLE SYMMETRIC DISTMAT WITHOUT CUSTOM LOOP OR USING SEQUENTIAL proxy LOOP
                dm <- quoted_call(
                    proxy::dist, x = x, y = NULL, method = distance, dots = dots
                )

                if (inherits(dm, "dist")) {
                    return(ret(dm))
                }
                else {
                    return(ret(base::as.matrix(dm), class = "crossdist"))
                }
            }
        }

        # WHOLE DISTMAT OR SUBDISTMAT OR NOT SYMMETRIC
        if (is.null(centroids)) centroids <- x
        dim_names <- list(names(x), names(centroids)) # x and centroids may change in parallel!
        x <- split_parallel(x)

        if (isTRUE(dots$pairwise)) {
            centroids <- split_parallel(centroids)
            validate_pairwise(x, centroids)
            combine <- c
        }
        else {
            centroids <- lapply(1L:foreach::getDoParWorkers(), function(dummy) { centroids })
            if (length(centroids) > length(x)) centroids <- centroids[1L:length(x)] # nocov
            combine <- rbind
        }

        d <- foreach(
            x = x, centroids = centroids,
            .combine = combine,
            .multicombine = TRUE,
            .packages = control$packages,
            .export = export
        ) %op% {
            if (!check_consistency(dist_entry$names[1L], "dist")) {
                do_call(proxy::pr_DB$set_entry, dist_entry)
            }

            quoted_call(proxy::dist, x = x, y = centroids, method = distance, dots = dots)
        }

        if (isTRUE(dots$pairwise)) {
            attr(d, "class") <- "pairdist"
        }
        else {
            attr(d, "class") <- "crossdist"
            attr(d, "dimnames") <- dim_names
        }
        # return
        ret(d)
    }

    # return closure
    distfun
}
