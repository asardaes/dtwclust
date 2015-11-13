# ========================================================================================================
# Return a custom family for kcca
# ========================================================================================================

kccaFamilies <- function(distance, cent, window.size, norm, distmat, ...) {

     ## Closures will capture the values of the constants,
     ## so that I don't have to use environments anymore...

     dots <- list(...)

     dtwdist <- switch(EXPR = distance,

                       ## Full DTW with L1 norm
                       dtw = {
                            foo <- function(x, centers) {
                                 x <- consistency_check(x, "tsmat")
                                 centers <- consistency_check(centers, "tsmat")

                                 if (!is.null(distmat))
                                      d <- dsub_pam(x, centers, distmat)
                                 else if (is.null(window.size)) {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "DTW", dist.method = "L1")
                                 } else {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "DTW", dist.method = "L1",
                                                       window.type = "slantedband",
                                                       window.size = window.size)
                                 }

                                 d
                            }

                            foo
                       },

                       ## Full DTW with L2 norm
                       dtw2 = {
                            foo <- function(x, centers) {
                                 x <- consistency_check(x, "tsmat")
                                 centers <- consistency_check(centers, "tsmat")

                                 if (!is.null(distmat))
                                      d <- dsub_pam(x, centers, distmat)
                                 else if (is.null(window.size)) {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "DTW2")
                                 } else {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "DTW2",
                                                       window.type = "slantedband",
                                                       window.size = window.size)
                                 }

                                 d
                            }

                            foo
                       },

                       ## DTW with aid of lower bounds
                       dtw_lb = {
                            window.size <- consistency_check(window.size, "window")

                            foo <- function(x, centers) {
                                 if (!is.null(distmat))
                                      d <- dsub_pam(x, centers, distmat)
                                 else
                                      d <- dtw_lb(x, centers, window.size, norm, error.check=FALSE)

                                 d
                            }

                            foo
                       },


                       ## Lemire's improved lower bound
                       lbi = {
                            window.size <- consistency_check(window.size, "window")

                            foo <- function(x, centers) {
                                 if (!is.null(distmat))
                                      d <- dsub_pam(x, centers, distmat)
                                 else {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "LBI",
                                                       norm = norm,
                                                       window.size = window.size,
                                                       error.check=FALSE)
                                 }

                                 d
                            }

                            foo
                       },

                       ## Keogh's lower bound
                       lbk = {
                            window.size <- consistency_check(window.size, "window")

                            foo <- function(x, centers) {
                                 if (!is.null(distmat))
                                      d <- dsub_pam(x, centers, distmat)
                                 else {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "LBK",
                                                       norm = norm,
                                                       window.size = window.size,
                                                       error.check=FALSE)
                                 }

                                 d
                            }

                            foo
                       },

                       ## Paparrizos' shape-based distance
                       sbd = {
                            foo <- function(x, centers) {
                                 x <- consistency_check(x, "tsmat")
                                 centers <- consistency_check(centers, "tsmat")

                                 if (!is.null(distmat))
                                      d <- dsub_pam(x, centers, distmat)
                                 else {
                                      d <- proxy::dist(x = x, y = centers,
                                                       method = "SBD", error.check = FALSE)
                                 }

                                 d
                            }

                            foo
                       },

                       ## Otherwise
                       foo <- function(x, centers) {
                            x <- consistency_check(x, "tsmat")
                            centers <- consistency_check(centers, "tsmat")

                            if (!is.null(distmat))
                                 d <- dsub_pam(x, centers, distmat)
                            else {
                                 ## do.call to ensure that the 'dots' argument is passed as '...'
                                 d <- do.call(proxy::dist,
                                              args = c(list(x = x,
                                                            y = centers,
                                                            method = distance),
                                                       dots))
                            }

                            d
                       }
     )

     family <- flexclust::kccaFamily(name = distance,

                                     dist = dtwdist,

                                     cent = cent)

     family
}

# ========================================================================================================
# Check consistency, used by other functions
# ========================================================================================================

consistency_check <- function(obj, case) {

     case <- match.arg(case, c("ts", "tslist", "vltslist", "window", "tsmat"))

     if (case == "ts") {
          if (!is.numeric(obj)) {
               stop("The series must be numeric")
          }
          if (!is.null(dim(obj)) && min(dim(obj)) != 1) {
               stop("The series must be univariate vectors")
          }
          if (length(obj) < 1) {
               stop("The series must have at least one point")
          }
          if (any(is.na(obj))) {
               stop("There are missing values in the series")
          }

     } else if (case == "tslist") {
          if (class(obj) != "list")
               stop("Series must be provided in the form of a list")

          if (length(obj) < 1)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) {
               ifelse(!is.null(dim(obj)) && min(dim(obj)) != 1, TRUE, FALSE)
          })))
               stop("Each element of the list must be a univariate vector")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (length(unique(sapply(obj, length))) > 1)
               stop("All series must have the same length")

          if (any(sapply(obj, length) < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "vltslist") {

          ## list of variable-length time series

          if (class(obj) != "list")
               stop("Series must be provided in the form of a list")

          if (length(obj) < 1)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) {
               ifelse(!is.null(dim(obj)) && min(dim(obj)) != 1, TRUE, FALSE)
          })))
               stop("Each element of the list must be a univariate vector")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (any(sapply(obj, length) < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "window") {
          if (is.null(obj)) {
               stop("Please provide the 'window.size' parameter")
          }
          if (obj <= 0) {
               stop("Window width must be larger than 0")
          }

          return(round(obj))

     } else if (case == "tsmat") {

          if (is.matrix(obj))
               obj <- lapply(seq_len(nrow(obj)), function(i) obj[i,])
          else if (is.numeric(obj))
               obj <- list(obj)
          else if (!is.list(obj))
               stop("Unsupported type for data")

          return(obj)

     } else {
          stop("Possibly a typo in consistency_check")
     }

     invisible(NULL)
}

# ========================================================================================================
# Custom functions to use for PAM centroids
# ========================================================================================================

# Custom function (closure) to obtain centroids when using PAM

allcent_pam <- function(distmat, distfun){
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
}

# Custom function to speed up the subsetting of the distance matrix for PAM

dsub_pam <- function(x, centers, distmat) {

     C <- consistency_check(centers, "tsmat")
     X <- consistency_check(x, "tsmat")

     indXC <- sapply(C, FUN = function(i.c) {
          i.row <- sapply(X, function(i.x) {
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

     d
}

# Pre-computing when using PAM. The whole distance matrix is calculated once and then reused

distmat_pam <- function(x, fam) {
     x <- consistency_check(x, "tsmat")

     d <- fam@dist(x, x)

     d
}

# ========================================================================================================
# Custom functions when using SBD averaging
# ========================================================================================================

# Custom function (closure) to obtain centroids when using SBD averaging

allcent_se <- function() {

     cluster.old <- NULL

     foo <- function(x, cluster, k) {

          # This will be read from parent environment
          cen <- get("centers", envir=parent.frame())
          C <- consistency_check(cen, "tsmat")

          if (is.matrix(x)) {
               X <- split.data.frame(x, cluster)
          } else if (is.list(x)) {
               X <- split(x, cluster)
               X <- lapply(X, function(l) t(sapply(l, rbind)))
          }

          cl <- sort(unique(cluster))

          ## Check which clusters changed
          if (is.null(cluster.old)) {
               indChanged <- cl

          } else {
               idchng <- cluster != cluster.old
               indChanged <- which( cl %in% sort(unique(c( cluster[idchng], cluster.old[idchng] ))) )
          }

          if(length(indChanged) == 0) {
               return(cen)
          }

          cluster.old <<- cluster # update

          new.C <- cen # initialize

          ## recompute centers for the clusters that changed
          if(is.matrix(cen)) {
               new.C[indChanged, ] <- t(mapply(X[indChanged], C[indChanged],
                                               FUN = function(x, c) {
                                                    new.c <- shape_extraction(x, c)

                                                    new.c
                                               }))

          } else if (is.list(cen)) {
               new.C[indChanged] <- mapply(X[indChanged], C[indChanged], SIMPLIFY = FALSE,
                                           FUN = function(x, c) {
                                                new.c <- shape_extraction(x, c)

                                                new.c
                                           })
          }

          new.C
     }

     foo
}

# Preprocessing when using shape_extraction (z-normalization)

preproc_se <- function (x) {
     if (is.matrix(x)) {
          xz <- apply(x, 1, zscore)
          t(xz)

     } else if (is.list(x)) {
          xz <- lapply(x, zscore)
     }
}

# ========================================================================================================
# Custom function to obtain centroids when using DBA
# ========================================================================================================

allcent_dba <- function(dba.iter, window.size, norm, trace) {

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
               indChanged <- which( cl %in% sort(unique(c( cluster[idchng], cluster.old[idchng] ))) )
          }

          if(length(indChanged) == 0) {
               return(cen)
          }

          cluster.old <<- cluster # update

          new.C <- cen # initialize

          ## Recompute centers for those clusters that changed
          if (is.matrix(cen)) {
               new.C[indChanged, ] <- t(mapply(X[indChanged], C[indChanged], indChanged,
                                               FUN = function(x, c, i) {
                                                    if(trace)
                                                         cat("\t\tComputing DBA center ", i, "...\n", sep = "")

                                                    new.c <- DBA(x, c,
                                                                 norm = norm, window.size = window.size,
                                                                 max.iter = dba.iter, error.check = FALSE)

                                                    new.c
                                               }))

          } else if (is.list(cen)) {
               new.C[indChanged] <- mapply(X[indChanged], C[indChanged], indChanged, SIMPLIFY = FALSE,
                                           FUN = function(x, c, i) {
                                                if(trace)
                                                     cat("\t\tComputing DBA center ", i, "...\n", sep = "")

                                                new.c <- DBA(x, c,
                                                             norm = norm, window.size = window.size,
                                                             max.iter = dba.iter, error.check = FALSE)

                                                new.c
                                           })
          }

          new.C
     }

     foo
}
