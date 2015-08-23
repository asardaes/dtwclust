# ========================================================================================================
# Return a custom family for kcca
# ========================================================================================================

kccaFamilies <- function(distance, cent, window.size, norm) {

     family <- switch(EXPR = distance,

                      ## Full DTW with L1 norm
                      dtw = flexclust::kccaFamily(name = "DTW",

                                                  dist = function(x, centers) {
                                                       window.size <- get("window.size", envir = attr(x, "env"))
                                                       distmat <- get("distmat", envir = attr(x, "env"))

                                                       x <- consistency_check(x, "tsmat")
                                                       centers <- consistency_check(centers, "tsmat")

                                                       if (!is.null(distmat))
                                                            d <- dsub_pam(x, centers)
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
                                                  },


                                                  cent = cent),

                      ## Full DTW with L2 norm
                      dtw2 = flexclust::kccaFamily(name = "DTW2",

                                                   dist = function(x, centers) {
                                                        window.size <- get("window.size", envir = attr(x, "env"))
                                                        distmat <- get("distmat", envir = attr(x, "env"))

                                                        x <- consistency_check(x, "tsmat")
                                                        centers <- consistency_check(centers, "tsmat")

                                                        if (!is.null(distmat))
                                                             d <- dsub_pam(x, centers)
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
                                                   },


                                                   cent = cent),


                      ## DTW with aid of lower bounds
                      dtw_lb = {
                           window.size <- consistency_check(window.size, "window")

                           flexclust::kccaFamily(name = "DTW_LB",

                                                 dist = function(x, centers) {
                                                      window.size <- get("window.size", envir = attr(x, "env"))
                                                      norm <- get("norm", envir = attr(x, "env"))
                                                      distmat <- get("distmat", envir = attr(x, "env"))

                                                      if (!is.null(distmat))
                                                           d <- dsub_pam(x, centers)
                                                      else
                                                           d <- dtw_lb(x, centers, window.size, norm, error.check=FALSE)

                                                      d
                                                 },

                                                 cent = cent)
                      },



                      ## Lemire's improved lower bound
                      lbi = {
                           window.size <- consistency_check(window.size, "window")

                           flexclust::kccaFamily(name = "LB_Improved",

                                                 dist = function(x, centers) {
                                                      window.size <- get("window.size", envir = attr(x, "env"))
                                                      norm <- get("norm", envir = attr(x, "env"))
                                                      distmat <- get("distmat", envir = attr(x, "env"))

                                                      if (!is.null(distmat))
                                                           d <- dsub_pam(x, centers)
                                                      else {
                                                           d <- proxy::dist(x = x, y = centers,
                                                                            method = "LBI",
                                                                            norm = norm,
                                                                            window.size = window.size,
                                                                            error.check=FALSE)
                                                      }

                                                      d
                                                 },


                                                 cent = cent)
                      },

                      ## Keogh's lower bound
                      lbk = {
                           window.size <- consistency_check(window.size, "window")

                           flexclust::kccaFamily(name = "LB_Keogh",

                                                 dist = function(x, centers) {
                                                      window.size <- get("window.size", envir = attr(x, "env"))
                                                      norm <- get("norm", envir = attr(x, "env"))
                                                      distmat <- get("distmat", envir = attr(x, "env"))

                                                      if (!is.null(distmat))
                                                           d <- dsub_pam(x, centers)
                                                      else {
                                                           d <- proxy::dist(x = x, y = centers,
                                                                            method = "LBK",
                                                                            norm = norm,
                                                                            window.size = window.size,
                                                                            error.check=FALSE)
                                                      }

                                                      d
                                                 },


                                                 cent = cent)
                      },

                      ## Paparrizos' shape-based distance
                      sbd = flexclust::kccaFamily(name = "SBD",

                                                  dist = function(x, centers) {
                                                       distmat <- get("distmat", envir = attr(x, "env"))

                                                       x <- consistency_check(x, "tsmat")
                                                       centers <- consistency_check(centers, "tsmat")

                                                       if (!is.null(distmat))
                                                            d <- dsub_pam(x, centers)
                                                       else {
                                                            d <- proxy::dist(x = x, y = centers,
                                                                             method = "SBD")
                                                       }

                                                       d
                                                  },


                                                  cent = cent),


                      ## Otherwise
                      flexclust::kccaFamily(name = distance,

                                            dist = function(x, centers) {
                                                 distmat <- get("distmat", envir = attr(x, "env"))

                                                 x <- consistency_check(x, "tsmat")
                                                 centers <- consistency_check(centers, "tsmat")

                                                 if (!is.null(distmat))
                                                      d <- dsub_pam(x, centers)
                                                 else {
                                                      d <- proxy::dist(x = x, y = centers,
                                                                       method = distance)
                                                 }

                                                 d
                                            },


                                            cent = cent)
     )

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
          if (!is.vector(obj)) {
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

          if (any(!sapply(obj, is.vector)))
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

          if (any(!sapply(obj, is.vector)))
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

          env <- attr(obj, "env")

          if (is.matrix(obj))
               obj <- lapply(seq_len(nrow(obj)), function(i) obj[i,])
          else if (is.numeric(obj))
               obj <- list(obj)
          else if (!is.list(obj))
               stop("Unsupported type for data")

          attr(obj, "env") <- env

          return(obj)

     } else {
          stop("Possibly a typo in consistency_check")
     }

     invisible(NULL)
}

# ========================================================================================================
# Custom functions to use for PAM centroids
# ========================================================================================================

# Custom function to obtain centroids when using PAM

allcent_pam <- function(x, cluster, k) {
     distmat <- get("distmat", envir = attr(x, "env"))

     indX <- lapply(sort(unique(cluster)), function(i) {
          which(cluster == i)
     })

     C <- sapply(indX, function(i.x) {
          d <- apply(distmat[i.x, i.x, drop=FALSE], 1, sum)

          if (is.matrix(x))
               x[i.x[which.min(d)], ]
          else if (is.list(x))
               x[i.x[which.min(d)]]
     })

     if (is.matrix(C))
          return(t(C))
     else
          return(C)
}

# Custom function to speed up the subsetting of the distance matrix for PAM

dsub_pam <- function(x, centers) {
     distmat <- get("distmat", envir = attr(x, "env"))

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

# Preprocessing when using PAM. The whole distance matrix is calculated once and reused

distmat_pam <- function(x, fam) {
     x <- consistency_check(x, "tsmat")

     d <- fam@dist(x, x)

     d
}

# ========================================================================================================
# Custom functions when using SBD averaging
# ========================================================================================================

# Custom function to obtain centroids when using SBD averaging

allcent_se <- function(x, cluster, k) {

     # This will be read from parent environment
     cen <- get("centers", envir=parent.frame())
     C <- consistency_check(cen, "tsmat")

     X <- split.data.frame(x, cluster)

     cl <- sort(unique(cluster))
     ncl <- length(cl)

     ## notice that the centers have to be in ascending order
     new.C <- mapply(X, C[cl],
                     FUN = function(x, c) {
                          new.c <- shape_extraction(x, c)

                          new.c
                     })

     t(new.C)
}

# Preprocessing when using shape_extraction (z-normalization)

preproc_se <- function (x) {
     xz <- apply(x, 1, zscore)
     attr(xz,"env") <- attr(x,"env")
     t(xz)
}

# ========================================================================================================
# Custom function to obtain centroids when using DBA
# ========================================================================================================

allcent_dba <- function(x, cluster, k) {
     max.iter <- get("dba.iter", envir = attr(x, "env"))

     # This will be read from parent environment
     cen <- get("centers", envir=parent.frame())
     C <- consistency_check(cen, "tsmat")

     cl <- sort(unique(cluster))

     if (is.matrix(x))
          X <- split.data.frame(x, cluster)
     else if (is.list(x)) {
          indXC <- lapply(cl, function(i.cl) which(cluster == i.cl))
          X <- lapply(indXC, function(ii) x[ii])
     } else
          stop("Invalid format for data")

     ## notice that the centers have to be in ascending order
     new.C <- mapply(X, C[cl],
                     FUN = function(x, c) {
                          new.c <- DBA(x, c, max.iter = max.iter, error.check = FALSE)

                          new.c
                     })

     if (is.matrix(new.C))
          return(t(new.C))
     else
          return(new.C)
}
