# ========================================================================================================
# Return a custom family for kcca
# ========================================================================================================

kccaFamilies <- function(distance, cent, window.size, norm, distmat, packages, ...) {

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

     dtwdist <- function(x, centers, whole = FALSE, ...) {

          ## Just in case, always empty for now
          dots2 <- c(dots, list(...))

          x <- consistency_check(x, "tsmat")
          centers <- consistency_check(centers, "tsmat")

          if (!is.null(distmat)) {
               ## distmat matrix already full, just subset it

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

          } else {
               ## Attempt to calculate in parallel?
               do_par <- pr_DB$get_entry(distance)$loop && foreach::getDoParRegistered()

               if (do_par) {

                    ## I need to re-register any custom distances in each parallel worker
                    dist_entry <- proxy::pr_DB$get_entry(distance)

                    ## variables from the parent environment that should be exported
                    export <- c("distance", "cent", "window.size", "norm", "distmat",
                                "dots", "window.type")

                    tasks <- foreach::getDoParWorkers()

                    if (whole) {
                         ## Whole distmat is calculated

                         ## by column so that I can assign it with upper.tri at the end
                         pairs <- call_pairs(length(x), byrow = FALSE)

                         tasks <- parallel::splitIndices(nrow(pairs), tasks)

                         pairs <- lapply(tasks, function(id){
                              pairs[id,]
                         })

                         d <- foreach(pairs = pairs,
                                      .combine = c,
                                      .multicombine = TRUE,
                                      .packages = c("dtwclust", packages),
                                      .export = export) %dopar% {

                                           if (!proxy::pr_DB$entry_exists(dist_entry$names[1]))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## Does the registered function posses '...' in its definition?
                                           has_dots <- is.function(dist_entry$FUN) &&
                                                any(grepl("...", names(as.list(args(dist_entry$FUN))),
                                                          fixed = TRUE))

                                           if (has_dots) {
                                                dd <- do.call(proxy::dist,
                                                              args = c(list(x = x[pairs[,1]],
                                                                            y = centers[pairs[,2]],
                                                                            method = distance,
                                                                            window.type = window.type,
                                                                            window.size = window.size,
                                                                            norm = norm,
                                                                            error.check = FALSE,
                                                                            pairwise = TRUE),

                                                                       dots2))
                                           } else {
                                                dd <- proxy::dist(x[pairs[,1]], centers[pairs[,2]],
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

                         x <- lapply(tasks, function(idx) {
                              x[idx]
                         })

                         d <- foreach(x = x,
                                      .combine = rbind,
                                      .multicombine = TRUE,
                                      .packages = c("dtwclust", packages),
                                      .export = export) %dopar% {

                                           if (!proxy::pr_DB$entry_exists(dist_entry$names[1]))
                                                do.call(proxy::pr_DB$set_entry, dist_entry)

                                           ## Does the registered function posses '...' in its definition?
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

                    ## Does the registered function posses '...' in its definition?
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

     family <- flexclust::kccaFamily(name = distance,

                                     dist = dtwdist,

                                     cent = cent)

     family
}

# ========================================================================================================
# Check consistency, used by other functions
# ========================================================================================================

#' Check consistency of provided object
#'
#' Used extensively internally. Exported just for detection in case of parallel computing.
#'
#' @param obj Object to check.
#' @param case One of \code{c("ts", "tslist", "vltslist", "window", "tsmat")}
#'
#' @return \code{NULL} invisibly
#'
#' @export

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

# Pre-computing when using PAM. The whole distance matrix is calculated once and then reused

distmat_pam <- function(x, fam) {
     x <- consistency_check(x, "tsmat")

     d <- fam@dist(x, x, TRUE)

     d
}

# Create combinations of all possible pairs

call_pairs <- function(n = 2L, byrow = TRUE) {
     if (n < 2)
          stop("I need at least two elements to create pairs.")

     .Call("pairs", n, byrow, PACKAGE = "dtwclust")
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
          xz <- zscore(x)
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
