# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

dtwdistfun <- function(distance, window.size, norm, distmat, packages, ...) {

     ## If I call this function is because 'distance' was a character

     dots <- list(...)

     needs_window <- c("dtw_lb", "lbk", "lbi")

     if (distance %in% needs_window)
          window.size <- consistency_check(window.size, "window")

     if (is.null(window.size))
          window.type <- "none"
     else
          window.type <- "slantedband"

     ## Closures will capture the values of the constants

     dtwdist <- function(x, centers = NULL, ...) {

          ## Just in case, always empty for now
          dots2 <- c(dots, list(...))

          x <- consistency_check(x, "tsmat")

          if (!is.null(centers))
               centers <- consistency_check(centers, "tsmat")

          if (!is.null(distmat)) {
               ## distmat matrix already calculated, just subset it

               indXC <- sapply(centers, FUN = function(i.c) {
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

               d <- distmat[ , indXC]

          } else {
               ## Attempt to calculate in parallel?
               do_par <- check_parallel(distance)

               if (do_par) {

                    ## I need to re-register any custom distances in each parallel worker
                    dist_entry <- proxy::pr_DB$get_entry(distance)

                    ## variables from the parent environment that should be exported
                    export <- c("distance", "window.size", "norm",
                                "window.type", "consistency_check")

                    tasks <- foreach::getDoParWorkers()

                    if (is.null(centers)) {
                         ## Whole distmat is calculated

                         ## by column so that I can assign it with upper.tri at the end
                         pairs <- call_pairs(length(x), byrow = FALSE)

                         tasks <- parallel::splitIndices(nrow(pairs), tasks)
                         tasks <- tasks[sapply(tasks, length, USE.NAMES = FALSE) != 0]

                         pairs <- lapply(tasks, function(id){
                              pairs[id,]
                         })

                         d <- foreach(pairs = pairs,
                                      .combine = c,
                                      .multicombine = TRUE,
                                      .packages = packages,
                                      .export = export) %dopar% {

                                           if (!proxy::pr_DB$entry_exists(dist_entry$names[1]))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## Does the registered function possess '...' in its definition?
                                           has_dots <- is.function(dist_entry$FUN) &&
                                                any(grepl("...", names(as.list(args(dist_entry$FUN))),
                                                          fixed = TRUE))

                                           if (has_dots) {
                                                dd <- do.call(proxy::dist,
                                                              args = c(list(x = x[pairs[,1]],
                                                                            y = x[pairs[,2]],
                                                                            method = distance,
                                                                            window.type = window.type,
                                                                            window.size = window.size,
                                                                            norm = norm,
                                                                            error.check = FALSE,
                                                                            pairwise = TRUE),

                                                                       dots2))
                                           } else {
                                                dd <- proxy::dist(x[pairs[,1]], x[pairs[,2]],
                                                                  method = distance, pairwise = TRUE)
                                           }

                                           dd

                                      }

                         D <- matrix(0, nrow = length(x), ncol = length(x))
                         D[upper.tri(D)] <- d
                         D <- t(D)
                         D[upper.tri(D)] <- d

                         d <- D
                         attr(d, "class") <- "crossdist"
                         attr(d, "method") <- distance

                    } else {
                         ## Only subset of distmat is calculated

                         tasks <- parallel::splitIndices(length(x), tasks)
                         tasks <- tasks[sapply(tasks, length, USE.NAMES = FALSE) != 0]

                         x <- lapply(tasks, function(idx) {
                              x[idx]
                         })

                         d <- foreach(x = x,
                                      .combine = rbind,
                                      .multicombine = TRUE,
                                      .packages = packages,
                                      .export = export) %dopar% {

                                           if (!proxy::pr_DB$entry_exists(dist_entry$names[1]))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## Does the registered function possess '...' in its definition?
                                           has_dots <- is.function(dist_entry$FUN) &&
                                                any(grepl("...", names(as.list(args(dist_entry$FUN))),
                                                          fixed = TRUE))

                                           if (has_dots) {
                                                dd <- do.call(proxy::dist,
                                                              args = c(list(x = x, y = centers,
                                                                            method = distance,
                                                                            window.type = window.type,
                                                                            window.size = window.size,
                                                                            norm = norm,
                                                                            error.check = FALSE),

                                                                       dots2))
                                           } else {
                                                dd <- proxy::dist(x, centers, method = distance)
                                           }

                                           dd

                                      }
                    }

               } else {
                    ## NO PARALLEL

                    if (is.null(centers))
                         centers <- x

                    ## Does the registered function possess '...' in its definition?
                    has_dots <- is.function(pr_DB$get_entry(distance)$FUN) &&
                         any(grepl("...", names(as.list(args(pr_DB$get_entry(distance)$FUN))),
                                   fixed = TRUE))

                    if (has_dots) {
                         ## If it has '...', put everything there and let it use whatever it needs

                         ## do.call to ensure that the 'dots2' argument is passed as '...'
                         d <- do.call(proxy::dist,
                                      args = c(list(x = x, y = centers,
                                                    method = distance,
                                                    window.type = window.type, window.size = window.size,
                                                    norm = norm, error.check = FALSE),

                                               dots2))
                    } else {
                         ## Otherwise just call it like this
                         d <- proxy::dist(x, centers, method = distance)
                    }

               }
          }

          d
     }

     dtwdist
}
