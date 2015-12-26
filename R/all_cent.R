# ========================================================================================================
# Custom functions to calculate centroids
# ========================================================================================================

all_cent <- function(case = NULL, distmat, distfun, control) {

     if (is.function(case))
          return(case)
     else if (!is.character(case))
          stop("Centroid definition must be either a function or a character")

     case <- match.arg(case, c("mean", "median", "shape", "dba", "pam"))

     pam_cent <- function(x, cluster, cl, Xsplit, ...) {

          if(is.null(distmat)) {
               C <- sapply(Xsplit, function(xsub) {
                    distmat <- distfun(xsub, xsub)

                    d <- apply(distmat, 1, sum)

                    xsub[which.min(d)]
               })

          } else {
               indX <- lapply(cl, function(id.cl) which(cluster == id.cl))

               C <- sapply(indX, function(i.x) {
                    d <- apply(distmat[i.x, i.x, drop=FALSE], 1L, sum)

                    x[i.x[which.min(d)]]
               })
          }

          C
     }

     shape_cent <- function(Xsplit, cluster, cen, cl, clustold, ...) {

          ## Check which clusters changed
          if (all(clustold == 0L)) {
               indChanged <- cl

          } else {
               idchng <- cluster != clustold
               indChanged <- union(cluster[idchng], clustold[idchng])
          }

          if(length(indChanged) == 0) {
               return(cen)
          }

          new.C <- cen # initialize

          ## recompute centers for the clusters that changed
          if (check_parallel()) {
               # in parallel
               indChanged <- split_parallel(indChanged, length(indChanged))

               exclude <- setdiff(ls(), c("Xsplit", "cen"))

               sub.C <- foreach(indChanged = indChanged,
                                .combine = c,
                                .multicombine = TRUE,
                                .packages = "dtwclust",
                                .noexport = exclude) %dopar% {
                                     mapply(Xsplit[indChanged], cen[indChanged],
                                            SIMPLIFY = FALSE,
                                            FUN = function(x, c) {
                                                 new.c <- shape_extraction(x, c)

                                                 new.c
                                            })
                                }
          } else {
               # not in parallel
               sub.C <- mapply(Xsplit[indChanged], cen[indChanged],
                               SIMPLIFY = FALSE,
                               FUN = function(x, c) {
                                    new.c <- shape_extraction(x, c)

                                    new.c
                               })
          }

          new.C[unlist(indChanged)] <- sub.C

          new.C
     }

     dba_cent <- function(Xsplit, cluster, cen, cl, clustold, ...) {

          ## Check which clusters changed
          if (all(clustold == 0)) {
               indChanged <- cl

          } else {
               idchng <- cluster != clustold
               indChanged <- union(cluster[idchng], clustold[idchng])
          }

          if(length(indChanged) == 0) {
               return(cen)
          }

          new.C <- cen # initialize

          ## Recompute centers for those clusters that changed
          if (check_parallel()) {
               # in parallel
               indChanged <- split_parallel(indChanged, length(indChanged))

               include <- c("control") # parent env
               exclude <- setdiff(ls(), c("Xsplit", "cen"))

               sub.C <- foreach(indChanged = indChanged,
                                .combine = c,
                                .multicombine = TRUE,
                                .packages = "dtwclust",
                                .export = include,
                                .noexport = exclude) %dopar% {
                                     mapply(Xsplit[indChanged],
                                            cen[indChanged],
                                            indChanged,
                                            SIMPLIFY = FALSE,
                                            FUN = function(x, c, i) {
                                                 if(control@trace)
                                                      cat("\t\tComputing DBA center ",
                                                          i,
                                                          "...\n",
                                                          sep = "")

                                                 new.c <- DBA(x, c,
                                                              norm = control@norm,
                                                              window.size = control@window.size,
                                                              max.iter = control@dba.iter,
                                                              error.check = FALSE)

                                                 new.c
                                            })
                                }
          } else {
               # not in parallel
               sub.C <- mapply(Xsplit[indChanged],
                               cen[indChanged],
                               indChanged,
                               SIMPLIFY = FALSE,
                               FUN = function(x, c, i) {
                                    if(control@trace)
                                         cat("\t\tComputing DBA center ",
                                             i,
                                             "...\n",
                                             sep = "")

                                    new.c <- DBA(x, c,
                                                 norm = control@norm,
                                                 window.size = control@window.size,
                                                 max.iter = control@dba.iter,
                                                 error.check = FALSE)

                                    new.c
                               })
          }

          new.C[unlist(indChanged)] <- sub.C

          new.C
     }

     mean_cent <- function(Xsplit, ...) {
          Xsplit <- lapply(Xsplit, function(xx) do.call(rbind, xx))

          centers <- lapply(Xsplit, colMeans)

          centers
     }

     median_cent <- function(Xsplit, ...) {
          Xsplit <- lapply(Xsplit, function(xx) do.call(rbind, xx))

          centers <- lapply(Xsplit, colMedians)

          centers
     }

     allcent <- function(x, cluster, k, cen, clustold, ...) {
          Xsplit <- split(x, cluster)

          cl <- sort(unique(cluster))

          centfun <- paste0(case, "_cent")

          cen <- do.call(centfun,
                         c(list(x = x, cluster = cluster,
                                k = k, cen = cen,
                                Xsplit = Xsplit, cl = cl,
                                clustold = clustold),
                           list(...)))

          ## Any empty clusters?
          empty_cluster <- k - length(cen)

          ## If so, initialize new clusters
          if (empty_cluster) {
               any_rep <- logical(empty_cluster)

               while(TRUE) {
                    extra_cen <- x[sample(length(x), empty_cluster)]

                    for (id_extra in 1L:empty_cluster) {
                         any_rep[id_extra] <- any(sapply(cen, function(i.center) {
                              if (length(i.center) != length(extra_cen[[id_extra]]))
                                   ret <- FALSE
                              else
                                   ret <- all(i.center == extra_cen[[id_extra]])

                              ret
                         }))
                    }

                    if (all(!any_rep))
                         break
               }

               cen <- c(cen, extra_cen)
          }

          cen
     }

     allcent
}
