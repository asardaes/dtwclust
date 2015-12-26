# ========================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ========================================================================================================

dtwdistfun <- function(distance, control, distmat) {

     needs_window <- c("dtw_lb", "lbk", "lbi")

     if (distance %in% needs_window)
          control@window.size <- consistency_check(control@window.size, "window")

     if (is.null(control@window.size))
          window.type <- "none"
     else
          window.type <- "slantedband"

     ## Closures will capture the values of the constants

     dtwdist <- function(x, centers = NULL, ...) {

          if (is.null(distmat)) {
               ## Extra distance parameters in case of parallel computation
               ## They can be for the function or for proxy::dist
               dots <- list(...)

               dots$window.type <- window.type
               dots$window.size <- control@window.size
               dots$norm <- control@norm
               dots$error.check <- FALSE

               ## I need to re-register any custom distances in each parallel worker
               dist_entry <- proxy::pr_DB$get_entry(distance)

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
          }

          if (!is.null(distmat)) {

               if (is.null(centers)) {
                    d <- distmat

               } else {
                    ## distmat matrix already calculated, just subset it
                    id_XC <- attr(centers, "id_cent")

                    if (is.null(id_XC) || any(is.na(id_XC)))
                         id_XC <- sapply(centers, FUN = function(i.c) {
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
               ## Attempt to calculate distmat in parallel?
               if (check_parallel(distance)) {

                    ## variables from the parent environment that should be exported
                    export <- c("distance", "control",
                                "window.type", "consistency_check")

                    if (is.null(centers)) {
                         ## Whole distmat is calculated

                         ## strict pairwise as in proxy::dist doesn't make sense here, but this case needs it
                         dots$pairwise <- TRUE

                         ## by column so that I can assign it with upper.tri at the end
                         pairs <- call_pairs(length(x), byrow = FALSE)

                         pairs <- split_parallel(pairs, nrow(pairs), 1L)

                         d <- foreach(pairs = pairs,
                                      .combine = c,
                                      .multicombine = TRUE,
                                      .packages = control@packages,
                                      .export = export) %dopar% {

                                           if (!consistency_check(dist_entry$names[1], "dist", silent = TRUE))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## 'dots' has all extra arguments that are valid
                                           dd <- do.call(proxy::dist,
                                                         c(dots,
                                                           list(x = x[pairs[,1]],
                                                                y = x[pairs[,2]],
                                                                method = distance)))

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

                         if (!is.null(dots$pairwise) && dots$pairwise) {
                              centers <- split_parallel(centers, length(centers))
                              combine <- c
                         } else {
                              centers <- lapply(1:length(x), function(dummy) centers)
                              combine <- rbind
                         }

                         d <- foreach(x = x, centers = centers,
                                      .combine = combine,
                                      .multicombine = TRUE,
                                      .packages = control@packages,
                                      .export = export) %dopar% {

                                           if (!consistency_check(dist_entry$names[1], "dist", silent = TRUE))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## 'dots' has all extra arguments that are valid
                                           dd <- do.call(proxy::dist,
                                                         c(dots,
                                                           list(x = x,
                                                                y = centers,
                                                                method = distance)))

                                           dd

                                      }
                    }

               } else {
                    ## NO PARALLEL

                    if (is.null(centers)) {
                         centers <- x
                         dots$pairwise <- FALSE # override possible TRUE
                    }

                    ## 'dots' has all extra arguments that are valid
                    d <- do.call(proxy::dist,
                                 c(dots,
                                   list(x = x,
                                        y = centers,
                                        method = distance)))
               }
          }

          d
     }

     dtwdist
}
