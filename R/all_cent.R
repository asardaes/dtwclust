# ========================================================================================================
# Custom functions to calculate centroids
# ========================================================================================================

all_cent <- function(case = NULL, distmat, distfun, control) {

     if (is.function(case))
          return(case)
     else if (!is.character(case))
          stop("Centroid definition must be either a function or a supported string")

     allcent <- switch(EXPR = case,

                       pam = {
                            foo <- function(x, cluster, k, cen, ...) {

                                 indX <- lapply(sort(unique(cluster)), function(i) {
                                      which(cluster == i)
                                 })

                                 if(is.null(distmat)) {
                                      C <- sapply(indX, function(i.x) {
                                           xsub <- x[i.x]
                                           distmat <- distfun(xsub, xsub)

                                           d <- apply(distmat, 1, sum)

                                           x[i.x[which.min(d)]]
                                      })

                                 } else {
                                      C <- sapply(indX, function(i.x) {
                                           d <- apply(distmat[i.x, i.x, drop=FALSE], 1, sum)

                                           x[i.x[which.min(d)]]
                                      })
                                 }

                                 C
                            }

                            foo
                       },

                       shape = {
                            cluster.old <- NULL

                            foo <- function(x, cluster, k, cen, ...) {

                                 X <- split(x, cluster)

                                 cl <- sort(unique(cluster))

                                 ## Check which clusters changed
                                 if (is.null(cluster.old)) {
                                      indChanged <- cl

                                 } else {
                                      idchng <- cluster != cluster.old
                                      indChanged <- which( cl %in%
                                                                sort(unique(c( cluster[idchng],
                                                                               cluster.old[idchng] ))) )
                                 }

                                 if(length(indChanged) == 0) {
                                      return(cen)
                                 }

                                 cluster.old <<- cluster # update

                                 new.C <- cen # initialize

                                 ## recompute centers for the clusters that changed
                                 if (check_parallel()) {
                                      # in parallel
                                      indChanged <- split_parallel(indChanged, length(indChanged))

                                      exclude <- setdiff(ls(), c("X", "cen"))

                                      sub.C <- foreach(indChanged = indChanged,
                                                       .combine = c,
                                                       .multicombine = TRUE,
                                                       .packages = "dtwclust",
                                                       .noexport = exclude) %dopar% {
                                                            mapply(X[indChanged], cen[indChanged],
                                                                   SIMPLIFY = FALSE,
                                                                   FUN = function(x, c) {
                                                                        new.c <- shape_extraction(x, c)

                                                                        new.c
                                                                   })
                                                       }
                                 } else {
                                      # not in parallel
                                      sub.C <- mapply(X[indChanged], cen[indChanged],
                                                      SIMPLIFY = FALSE,
                                                      FUN = function(x, c) {
                                                           new.c <- shape_extraction(x, c)

                                                           new.c
                                                      })
                                 }

                                 new.C[unlist(indChanged)] <- sub.C

                                 new.C
                            }

                            foo
                       },

                       dba = {
                            cluster.old <- NULL

                            ## Closure
                            foo <- function(x, cluster, k, cen, ...) {

                                 cl <- sort(unique(cluster))

                                 indXC <- lapply(cl, function(i.cl) which(cluster == i.cl))
                                 X <- lapply(indXC, function(ii) x[ii])

                                 ## Check which clusters changed
                                 if (is.null(cluster.old)) {
                                      indChanged <- cl

                                 } else {
                                      idchng <- cluster != cluster.old
                                      indChanged <- which( cl %in%
                                                                sort(unique(c( cluster[idchng],
                                                                               cluster.old[idchng] ))) )
                                 }

                                 if(length(indChanged) == 0) {
                                      return(cen)
                                 }

                                 cluster.old <<- cluster # update

                                 new.C <- cen # initialize

                                 ## Recompute centers for those clusters that changed
                                 if (check_parallel()) {
                                      # in parallel
                                      indChanged <- split_parallel(indChanged, length(indChanged))

                                      include <- c("control") # parent env
                                      exclude <- setdiff(ls(), c("X", "cen"))

                                      sub.C <- foreach(indChanged = indChanged,
                                                       .combine = c,
                                                       .multicombine = TRUE,
                                                       .packages = "dtwclust",
                                                       .export = include,
                                                       .noexport = exclude) %dopar% {
                                                            mapply(X[indChanged],
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
                                      sub.C <- mapply(X[indChanged],
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

                            foo
                       },

                       mean = {
                            foo <- function(x, cluster, k, cen, ...) {
                                 X <- split(x, cluster)
                                 X <- lapply(X, function(xx) do.call(rbind, xx))

                                 centers <- lapply(X, colMeans)

                                 centers
                            }

                            foo
                       },

                       median = {
                            colMedians <- function(mat) { apply(mat, 2, stats::median) }

                            foo <- function(x, cluster, k, cen, ...) {
                                 X <- split(x, cluster)
                                 X <- lapply(X, function(xx) do.call(rbind, xx))

                                 centers <- lapply(X, colMedians)

                                 centers
                            }

                            foo
                       },

                       ## Otherwise
                       stop("Invalid centroid definition"))

     allcent
}
