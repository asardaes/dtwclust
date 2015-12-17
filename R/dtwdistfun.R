# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

dtwdistfun <- function(distance, control, distmat, ...) {

     ## If I call this function is because 'distance' was a character

     dots <- list(...)

     needs_window <- c("dtw_lb", "lbk", "lbi")

     if (distance %in% needs_window)
          control@window.size <- consistency_check(control@window.size, "window")

     if (is.null(control@window.size))
          window.type <- "none"
     else
          window.type <- "slantedband"

     ## Closures will capture the values of the constants

     dtwdist <- function(x, centers = NULL, ...) {

          ## Just in case, always empty for now
          dots2 <- c(dots, list(...))

          if (!is.null(distmat)) {

               if (is.null(centers)) {
                    d <- distmat

               } else {
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

                    d <- distmat[ , indXC, drop = FALSE]
               }

          } else {
               ## Attempt to calculate distmat in parallel?
               do_par <- check_parallel(distance)

               if (do_par) {

                    ## I need to re-register any custom distances in each parallel worker
                    dist_entry <- proxy::pr_DB$get_entry(distance)

                    ## variables from the parent environment that should be exported
                    export <- c("distance", "control",
                                "window.type", "consistency_check")

                    if (is.null(centers)) {
                         ## Whole distmat is calculated

                         ## by column so that I can assign it with upper.tri at the end
                         pairs <- call_pairs(length(x), byrow = FALSE)

                         pairs <- split_parallel(pairs, nrow(pairs), 1L)

                         d <- foreach(pairs = pairs,
                                      .combine = c,
                                      .multicombine = TRUE,
                                      .packages = control@packages,
                                      .export = export) %dopar% {

                                           if (!proxy::pr_DB$entry_exists(dist_entry$names[1]))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## Does the registered function possess '...' in its definition?
                                           has_dots <- is.function(dist_entry$FUN) && !is.null(formals(dist_entry$FUN)$...)

                                           if (has_dots) {
                                                dd <- do.call(proxy::dist,
                                                              args = c(list(x = x[pairs[,1]],
                                                                            y = x[pairs[,2]],
                                                                            method = distance,
                                                                            window.type = window.type,
                                                                            window.size = control@window.size,
                                                                            norm = control@norm,
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
                         x <- split_parallel(x, length(x))

                         d <- foreach(x = x,
                                      .combine = rbind,
                                      .multicombine = TRUE,
                                      .packages = control@packages,
                                      .export = export) %dopar% {

                                           if (!proxy::pr_DB$entry_exists(dist_entry$names[1]))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## Does the registered function possess '...' in its definition?
                                           has_dots <- is.function(dist_entry$FUN) && !is.null(formals(dist_entry$FUN)$...)

                                           if (has_dots) {
                                                dd <- do.call(proxy::dist,
                                                              args = c(list(x = x, y = centers,
                                                                            method = distance,
                                                                            window.type = window.type,
                                                                            window.size = control@window.size,
                                                                            norm = control@norm,
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
                    has_dots <- is.function(pr_DB$get_entry(distance)$FUN) && !is.null(formals(pr_DB$get_entry(distance)$FUN)$...)

                    if (has_dots) {
                         ## If it has '...', put everything there and let it use whatever it needs

                         ## do.call to ensure that the 'dots2' argument is passed as '...'
                         d <- do.call(proxy::dist,
                                      args = c(list(x = x, y = centers,
                                                    method = distance,
                                                    window.type = window.type, window.size = control@window.size,
                                                    norm = control@norm, error.check = FALSE),

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
