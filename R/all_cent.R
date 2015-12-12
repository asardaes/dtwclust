# ========================================================================================================
# Custom functions to calculate centroids
# ========================================================================================================

all_cent <- function(case = NULL,
                     distmat, distfun,
                     dba.iter, window.size, norm, trace) {

     if (is.function(case))
          return(case)
     else if (!is.character(case))
          stop("Possible mistake in all_cent.R")

     allcent <- switch(EXPR = case,

                       pam = {
                            foo <- function(x, cluster, k) {

                                 indX <- lapply(sort(unique(cluster)), function(i) {
                                      which(cluster == i)
                                 })

                                 if(is.null(distmat)) {
                                      C <- sapply(indX, function(i.x) {
                                           if (is.matrix(x)) {
                                                xsub <- x[i.x, , drop = FALSE] # function will deal with matrices
                                                distmat <- distfun(xsub, xsub)

                                           } else if (is.list(x)) {
                                                xsub <- x[i.x]
                                                distmat <- distfun(xsub, xsub)
                                           }

                                           d <- apply(distmat, 1, sum)

                                           if (is.matrix(x))
                                                x[i.x[which.min(d)], ]
                                           else if (is.list(x))
                                                x[i.x[which.min(d)]]
                                      })

                                 } else {
                                      C <- sapply(indX, function(i.x) {
                                           d <- apply(distmat[i.x, i.x, drop=FALSE], 1, sum)

                                           if (is.matrix(x))
                                                x[i.x[which.min(d)], ]
                                           else if (is.list(x))
                                                x[i.x[which.min(d)]]
                                      })
                                 }

                                 if (is.matrix(C))
                                      return(t(C))
                                 else
                                      return(C)
                            }

                            foo
                       },

                       shape = {
                            cluster.old <- NULL

                            foo <- function(x, cluster, k) {

                                 # This will be read from parent environment
                                 cen <- get("centers", envir=parent.frame())
                                 C <- consistency_check(cen, "tsmat")

                                 if (is.matrix(x))
                                      X <- split.data.frame(x, cluster)
                                 else if (is.list(x))
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
                                 if (foreach::getDoParRegistered() && foreach::getDoParWorkers() > 1L) {
                                      # in parallel
                                      tasks <- parallel::splitIndices(length(indChanged),
                                                                      foreach::getDoParWorkers())
                                      tasks <- tasks[sapply(tasks, length, USE.NAMES = FALSE) != 0]

                                      indChanged <- lapply(tasks, function(id) {
                                           indChanged[id]
                                      })

                                      exclude <- setdiff(ls(), c("X", "C", "cen"))

                                      sub.C <- foreach(indChanged = indChanged,
                                                       .combine = ifelse(is.matrix(cen), cbind, c),
                                                       .multicombine = TRUE,
                                                       .packages = "dtwclust",
                                                       .noexport = exclude) %dopar% {
                                                            mapply(X[indChanged], C[indChanged],
                                                                   SIMPLIFY = is.matrix(cen),
                                                                   FUN = function(x, c) {
                                                                        new.c <- shape_extraction(x, c)

                                                                        new.c
                                                                   })
                                                       }
                                 } else {
                                      # not in parallel
                                      sub.C <- mapply(X[indChanged], C[indChanged],
                                                      SIMPLIFY = is.matrix(cen),
                                                      FUN = function(x, c) {
                                                           new.c <- shape_extraction(x, c)

                                                           new.c
                                                      })
                                 }

                                 if (is.matrix(cen))
                                      new.C[unlist(indChanged), ] <- t(sub.C)
                                 else
                                      new.C[unlist(indChanged)] <- sub.C

                                 new.C
                            }

                            foo
                       },

                       dba = {
                            cluster.old <- NULL

                            ## Closure
                            foo <- function(x, cluster, k) {

                                 # This will be read from parent environment
                                 cen <- get("centers", envir=parent.frame())
                                 C <- consistency_check(cen, "tsmat")
                                 x <- consistency_check(x, "tsmat")

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
                                 if (foreach::getDoParRegistered() && foreach::getDoParWorkers() > 1L) {
                                      # in parallel
                                      tasks <- parallel::splitIndices(length(indChanged),
                                                                      foreach::getDoParWorkers())
                                      tasks <- tasks[sapply(tasks, length, USE.NAMES = FALSE) != 0]

                                      indChanged <- lapply(tasks, function(id) {
                                           indChanged[id]
                                      })

                                      include <- c("trace", "norm", "window.size", "dba.iter") # parent env
                                      exclude <- setdiff(ls(), c("X", "C", "cen"))

                                      sub.C <- foreach(indChanged = indChanged,
                                                       .combine = ifelse(is.matrix(cen), cbind, c),
                                                       .multicombine = TRUE,
                                                       .packages = "dtwclust",
                                                       .export = include,
                                                       .noexport = exclude) %dopar% {
                                                            mapply(X[indChanged],
                                                                   C[indChanged],
                                                                   indChanged,
                                                                   SIMPLIFY = is.matrix(cen),
                                                                   FUN = function(x, c, i) {
                                                                        if(trace)
                                                                             cat("\t\tComputing DBA center ",
                                                                                 i,
                                                                                 "...\n",
                                                                                 sep = "")

                                                                        new.c <- DBA(x, c,
                                                                                     norm = norm,
                                                                                     window.size = window.size,
                                                                                     max.iter = dba.iter,
                                                                                     error.check = FALSE)

                                                                        new.c
                                                                   })
                                                       }
                                 } else {
                                      # not in parallel
                                      sub.C <- mapply(X[indChanged],
                                                      C[indChanged],
                                                      indChanged,
                                                      SIMPLIFY = is.matrix(cen),
                                                      FUN = function(x, c, i) {
                                                           if(trace)
                                                                cat("\t\tComputing DBA center ",
                                                                    i,
                                                                    "...\n",
                                                                    sep = "")

                                                           new.c <- DBA(x, c,
                                                                        norm = norm,
                                                                        window.size = window.size,
                                                                        max.iter = dba.iter,
                                                                        error.check = FALSE)

                                                           new.c
                                                      })
                                 }

                                 if (is.matrix(cen))
                                      new.C[unlist(indChanged), ] <- t(sub.C)
                                 else
                                      new.C[unlist(indChanged)] <- sub.C

                                 new.C
                            }

                            foo
                       },

                       mean = {
                            foo <- function(x, cluster, k) {
                                 if (is.list(x)) {
                                      X <- split(x, cluster)
                                      X <- lapply(X, function(xx) do.call(rbind, xx))

                                      centers <- lapply(X, foo)

                                 } else {
                                      X <- split.data.frame(x, cluster)

                                      centers <- lapply(X, colMeans)
                                      centers <- do.call(rbind, centers)
                                 }

                                 centers
                            }

                            foo
                       },

                       median = {
                            colMedians <- function(mat) { apply(mat, 2, stats::median) }

                            foo <- function(x, cluster, k) {
                                 if (is.list(x)) {
                                      X <- split(x, cluster)
                                      X <- lapply(X, function(xx) do.call(rbind, xx))

                                      centers <- lapply(X, foo)

                                 } else {
                                      X <- split.data.frame(x, cluster)

                                      centers <- lapply(X, colMedians)
                                      centers <- do.call(rbind, centers)
                                 }

                                 centers
                            }

                            foo
                       },

                       ## Otherwise
                       stop("Invalid centroid definition"))

     allcent
}
