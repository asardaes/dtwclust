# ========================================================================================================
# Custom functions to calculate centroids
# ========================================================================================================

all_cent <- function(case = NULL, distmat, distfun, control, fuzzy = FALSE) {

     if (is.function(case))
          return(case)
     else if (!is.character(case))
          stop("Centroid definition must be either a function or a character")

     pam_cent <- function(x, x_split, cl_id, id_changed, ...) {

          if(is.null(distmat)) {
               new_cent <- lapply(x_split, function(xsub) {
                    distmat <- distfun(xsub, xsub)

                    d <- apply(distmat, 1L, sum)

                    i_cent <- xsub[[which.min(d)]]
                    attr(i_cent, "id_cent") <- pmatch(names(xsub[which.min(d)]), names(x))

                    i_cent
               })

          } else {
               id_x <- lapply(id_changed, function(cl_num) which(cl_id == cl_num))

               new_cent <- lapply(id_x, function(i_x) {
                    d <- apply(distmat[i_x, i_x, drop=FALSE], 1L, sum)

                    i_cent <- x[[i_x[which.min(d)]]]
                    attr(i_cent, "id_cent") <- i_x[which.min(d)]

                    i_cent
               })
          }

          new_cent
     }

     shape_cent <- function(x_split, cent, id_changed, cl_id, cl_old, ...) {

          check_parallel()

          x_split <- split_parallel(x_split)
          cent <- split_parallel(cent)

          new_cent <- foreach(x_split = x_split,
                              cent = cent,
                              .combine = c,
                              .multicombine = TRUE,
                              .packages = "dtwclust") %dopar% {
                                   mapply(x_split, cent,
                                          SIMPLIFY = FALSE,
                                          FUN = function(x, c) {
                                               new_c <- shape_extraction(x, c)

                                               new_c
                                          })
                              }

          new_cent
     }

     dba_cent <- function(x_split, cent, id_changed, cl_id, cl_old, ...) {

          check_parallel()

          x_split <- split_parallel(x_split)
          cent <- split_parallel(cent)

          new_cent <- foreach(x_split = x_split,
                              cent = cent,
                              .combine = c,
                              .multicombine = TRUE,
                              .packages = "dtwclust",
                              .export = "control") %dopar% {
                                   mapply(x_split, cent,
                                          SIMPLIFY = FALSE,
                                          FUN = function(x, c) {
                                               new_c <- DBA(x, c,
                                                            norm = control@norm,
                                                            window.size = control@window.size,
                                                            max.iter = control@dba.iter,
                                                            delta = control@delta,
                                                            error.check = FALSE)

                                               new_c
                                          })
                              }

          new_cent
     }

     mean_cent <- function(x_split, ...) {
          x_split <- lapply(x_split, function(xx) do.call(rbind, xx))

          new_cent <- lapply(x_split, colMeans)

          new_cent
     }

     median_cent <- function(x_split, ...) {
          x_split <- lapply(x_split, function(xx) do.call(rbind, xx))

          new_cent <- lapply(x_split, colMedians)

          new_cent
     }

     if (fuzzy) {
          allcent <- function(x, cl_id, k, cent, cl_old, ...) {
               u <- cl_id ^ control@fuzziness

               cent <- t(u) %*% do.call(rbind, x)
               cent <- apply(cent, 2L, "/", e2 = colSums(u))

               # Coerce back to list
               consistency_check(cent, "tsmat")
          }
     } else {
          allcent <- function(x, cl_id, k, cent, cl_old, ...) {
               ## Check which clusters changed
               if (all(cl_old == 0L)) {
                    id_changed <- sort(unique(cl_id))

               } else {
                    id_changed <- cl_id != cl_old
                    id_changed <- union(cl_id[id_changed], cl_old[id_changed])
               }

               if(length(id_changed) == 0L) {
                    return(cent)
               }

               ## Split data according to cluster memebership
               x_split <- split(x, factor(cl_id, levels = 1L:k))

               ## In case of empty new clusters
               empty_clusters <- which(lengths(x_split) == 0L)
               id_changed <- setdiff(id_changed, empty_clusters)

               ## Calculate new centers
               new_cent <- do.call(paste0(case, "_cent"),
                                   c(list(x = x,
                                          x_split = x_split[id_changed],
                                          cent = cent[id_changed],
                                          id_changed = id_changed,
                                          cl_id = cl_id,
                                          cl_old = cl_old),
                                     list(...)))

               cent[id_changed] <- new_cent

               ## Any empty clusters?
               num_empty <- length(empty_clusters)

               ## If so, initialize new clusters
               if (num_empty > 0L) {
                    ## Make sure no center is repeated (especially in case of PAM)
                    any_rep <- logical(num_empty)

                    while(TRUE) {
                         id_cent_extra <- sample(length(x), num_empty)
                         extra_cent <- x[id_cent_extra]

                         for (id_extra in 1L:num_empty) {
                              any_rep[id_extra] <- any(sapply(cent, function(i.center) {
                                   if (length(i.center) != length(extra_cent[[id_extra]]))
                                        ret <- FALSE
                                   else
                                        ret <- all(i.center == extra_cent[[id_extra]])

                                   ret
                              }))

                              if (case == "pam")
                                   attr(extra_cent[[id_extra]], "id_cent") <- id_cent_extra[id_extra]
                         }

                         if (all(!any_rep))
                              break
                    }

                    cent[empty_clusters] <- extra_cent
               }

               if (case == "pam")
                    attr(cent, "id_cent") <- sapply(cent, attr, which = "id_cent")

               cent
          }
     }

     allcent
}
