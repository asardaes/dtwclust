# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

ddist <- function(distance, control, distmat) {
     needs_window <- c("dtw_lb", "lbk", "lbi")

     if (distance %in% needs_window)
          control@window.size <- consistency_check(control@window.size, "window")

     if (is.null(control@window.size))
          window.type <- "none"
     else
          window.type <- "slantedband"

     ## Closures will capture the values of the constants

     distfun <- function(x, centroids = NULL, ...) {
          if (!is.null(list(...)$centers)) {
               warning("The 'centers' argument has been deprecated, please use 'centroids' instead.")

               if (is.null(centroids)) centroids <- list(...)$centers
          }

          if (!is.null(distmat)) {
               if (is.null(centroids)) {
                    d <- distmat

               } else {
                    ## distmat matrix already calculated, just subset it
                    id_XC <- attr(centroids, "id_cent")

                    if (is.null(id_XC) || any(is.na(id_XC)))
                         id_XC <- sapply(centroids, FUN = function(i.c) {
                              i.row <- sapply(x, function(i.x) {
                                   if (length(i.x) == length(i.c))
                                        ret <- all(i.x == i.c)
                                   else
                                        ret <- FALSE

                                   ret
                              })

                              ## Take the first one in case a series is repeated more than once in the dataset
                              which(i.row)[1]
                         })

                    d <- distmat[ , id_XC, drop = FALSE]
               }

          } else {
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

               ## Does the registered function possess '...' in its definition?
               has_dots <- is.function(dist_entry$FUN) && !is.null(formals(dist_entry$FUN)$...)

               ## If it doesn't, remove invalid arguments from 'dots'
               if (!has_dots) {
                    if (is.function(dist_entry$FUN))
                         valid_args <- union(names(formals(proxy::dist)), names(formals(dist_entry$FUN)))
                    else
                         valid_args <- names(formals(proxy::dist))

                    dots <- dots[intersect(names(dots), valid_args)]
               }

               ## Register doSEQ if necessary
               check_parallel()

               ## variables/functions from the parent environment that should be exported
               export <- c("distance", "consistency_check", "enlist")

               if (is.null(centroids) && control@symmetric && dist_entry$loop) {
                    ## WHOLE SYMMETRIC DISTMAT
                    ## Only half of it is computed

                    ## strict pairwise as in proxy::dist doesn't make sense here, but this case needs it
                    dots$pairwise <- TRUE

                    pairs <- call_pairs(length(x), lower = FALSE)

                    pairs <- split_parallel(pairs, 1L)

                    d <- foreach(pairs = pairs,
                                 .combine = c,
                                 .multicombine = TRUE,
                                 .packages = control@packages,
                                 .export = export) %dopar% {

                                      if (!consistency_check(dist_entry$names[1], "dist"))
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

               } else {
                    ## WHOLE OR SUBDISTMAT, NOT SYMMETRIC OR loop = FALSE

                    if (is.null(centroids)) {
                         centroids <- x

                         ## for dtw_basic_proxy
                         if (!check_parallel() && has_dots && (is.null(control@window.size) || !check_lengths(x)))
                              dots$symmetric <- TRUE
                    }

                    dim_names <- list(names(x), names(centroids))

                    if (!is.null(dots$pairwise) && dots$pairwise) {
                         if (length(x) != length(centroids))
                              stop("Both sets of data must have the same amount of series for pairwise calculation")

                         centroids <- split_parallel(centroids)
                         combine <- c

                    } else {
                         centroids <- lapply(1L:foreach::getDoParWorkers(), function(dummy) centroids)
                         combine <- rbind
                    }

                    x <- split_parallel(x)

                    d <- foreach(x = x, centroids = centroids,
                                 .combine = combine,
                                 .multicombine = TRUE,
                                 .packages = control@packages,
                                 .export = export) %dopar% {

                                      if (!consistency_check(dist_entry$names[1L], "dist"))
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
