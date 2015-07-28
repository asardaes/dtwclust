# ========================================================================================================
# Check consistency, used by other functions
# ========================================================================================================

consistency_check <- function(obj, case) {
     case <- match.arg(case,
                       c("ts", "tslist"))

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
          ## For now, only support for list of time series
          if (class(obj) != "list")
               stop("Series must be provided in the form of a list")

          if (length(obj) < 1)
               stop("List is empty")

          if (any(!sapply(obj, is.vector)))
               stop("Each element of the list must be a univariate vector")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (any(length(obj[[1]]) != sapply(obj, length)))
               stop("All series must have the same length")

          if (any(sapply(obj, length) < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")
     }
}

# ========================================================================================================
# Custom functions to use for PAM centroids
# ========================================================================================================

# Custom function to obtain centroids when using PAM

allcent_pam <- function(x, cluster, k) {
     distmat <- attr(x, "distmat")

     indX <- lapply(sort(unique(cluster)), function(i) {
          which(cluster %in% i)
     })

     C <- sapply(indX, function(i.x) {
          d <- apply(distmat[i.x, i.x, drop=FALSE], 1, sum)

          return( x[i.x[which.min(d)], ] )
     })

     return(t(C))
}

# Custom function to speed up the subsetting of the distance matrix for PAM

dsub_pam <- function(x, centers) {
     distmat <- attr(x, "distmat")

     C <- lapply(seq_len(nrow(centers)), function(i) centers[i,])

     indXC <- sapply(C, FUN = function(i.c) {
          i.row <- apply(x, 1, function(i.x) {
               all(i.x == i.c)
          })

          return(which(i.row)[1])
     })

     d <- distmat[ , indXC]

     return(d)
}

# Preprocessing when using PAM. The whole distance matrix is calculated once and assigned as an attribute
# of the data.

preproc_pam <- function(x) {
     fam <- get("family", envir=parent.frame())

     d <- fam@dist(x, x)

     attr(x, "distmat") <- d

     return(x)
}

# ========================================================================================================
# Custom functions when using SBD averaging
# ========================================================================================================

# Custom function to obtain centroids when using SBD averaging

allcent_se <- function(x, cluster, k) {
     # This will be read from parent environment
     cen <- get("centers", envir=parent.frame())
     C <- lapply(seq_len(nrow(cen)), function(i) cen[i,])

     X <- split.data.frame(x, cluster)

     cl <- sort(unique(cluster))
     ncl <- length(cl)

     new.C <- mapply(X, C[cl],
                     FUN = function(x, c) {
                          new.c <- shape_extraction(x, c)

                          return(new.c)
                     })

     new.C <- t(new.C)

     ## I'm not sure if this can happen
     if (ncl < k) {
          miss.ind <- !((1:k) %in% cl)
          miss.ind <- (1:k)[miss.ind]

          #new.C <- rbind(new.C, cen[miss.ind, ])
          new.C <- rbind( new.C, matrix(0, length(miss.ind), ncol(new.C)) )

          order.ind <- sort(c(cl, miss.ind), index.return=T)
          new.C <- new.C[order.ind$ix, ]
     }

     return(new.C)
}

# Preprocessing when using SBD averaging (z-normalization)

preproc_se <- function (x) {
     xz <- apply(x, 1, zscore)
     return(t(xz))
}

# ========================================================================================================
# Custom function to obtain centroids when using DBA
# ========================================================================================================

allcent_dba <- function(x, cluster, k) {
     # This will be read from parent environment
     cen <- get("centers", envir=parent.frame())
     C <- lapply(seq_len(nrow(cen)), function(i) cen[i,])

     X <- split.data.frame(x, cluster)

     cl <- sort(unique(cluster))
     #ncl <- length(cl)

     new.C <- mapply(X, C[cl],
                     FUN = function(x, c) {
                          new.c <- DBA(x, c, error.check = FALSE)

                          new.c
                     })

     new.C <- t(new.C)

     new.C
}
