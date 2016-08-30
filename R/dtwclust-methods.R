#' @title Methods for \code{dtwclust}
#'
#' @description
#' Methods associated with \code{\link{dtwclust-class}} objects.
#'
#' @details
#' Supported generics from the \code{flexclust} package are: \code{\link[flexclust]{randIndex}} and
#' \code{\link[flexclust]{clusterSim}}.
#'
#' @name dtwclust-methods
#' @rdname dtwclust-methods
#' @include dtwclust-classes.R
#'
#' @seealso \code{\link{dtwclust-class}}, \code{\link{dtwclust}}, \code{\link[ggplot2]{ggplot}}
#'
NULL

# ========================================================================================================
# Custom initialize to avoid infinite recursion
# See https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16629
# ========================================================================================================

setMethod("initialize", "dtwclust",
          function(.Object, ..., call) {
               .Object <- callNextMethod(.Object = .Object, ...)

               if(!missing(call))
                    .Object@call <- call

               .Object
          })

# ========================================================================================================
# Show
# ========================================================================================================

#' @details
#' Show method displays basic information from the clustering results.
#'
#' @rdname dtwclust-methods
#' @aliases show,dtwclust-method
#'
#' @param object,x An object of class \code{\link{dtwclust-class}} as returned by \code{\link{dtwclust}}.
#'
#' @exportMethod show
#'
setMethod("show", "dtwclust",
          function(object) {
               print(object@call)
               cat("\n")

               cat(object@type, "clustering with", object@k, "clusters\n")
               cat("Using", object@distance, "distance\n")

               if (object@type == "partitional")
                    cat("Using", object@centroid, "centroids\n")
               if (object@type == "hierarchical")
                    cat("Using", object@method, "linkage\n")
               if (object@preproc != "none")
                    cat("Using", object@preproc, "preprocessing\n")

               cat("\nTime required for analysis:\n")
               print(object@proctime)

               if (object@type == "fuzzy") {
                    cat("\nHead of fuzzy memberships:\n\n")
                    print(utils::head(object@fcluster))

               } else {
                    cat("\nCluster sizes with average intra-cluster distance:\n\n")
                    print(object@clusinfo)
               }
          })

# ========================================================================================================
# update from stats
# ========================================================================================================

#' @details
#' The \code{update} method takes the original function call, replaces any provided argument and optionally
#' evaluates the call again. Use \code{evaluate = FALSE} if you want to get the
#' unevaluated call.
#'
#' @rdname dtwclust-methods
#' @aliases update,dtwclust-method
#'
#' @param evaluate Logical. Defaults to \code{TRUE} and evaluates the updated call, which will result in
#' a new \code{dtwclust} object. Otherwise, it returns the unevaluated call.
#'
#' @exportMethod update
#'
setMethod("update", "dtwclust",
          function(object, ..., evaluate = TRUE) {

               args <- as.pairlist(list(...))

               if (length(args) == 0L) {
                    message("Nothing to be updated")

                    if (evaluate)
                         return(object)
                    else
                         return(object@call)
               }

               new_call <- object@call
               new_call[names(args)] <- args

               if (evaluate)
                    ret <- eval.parent(new_call, n = 2L)
               else
                    ret <- new_call

               ret
          })

# ========================================================================================================
# predict from stats
# ========================================================================================================

#' @details
#' The \code{predict} generic can take the usual \code{newdata} argument and it returns the cluster(s) to which
#' the data belongs; if \code{NULL}, it simply returns the obtained cluster indices. It preprocesses
#' the data with the corresponding function if available.
#'
#' @rdname dtwclust-methods
#' @aliases predict,dtwclust-method
#'
#' @param newdata New data to be evaluated. It can take any of the supported formats of \code{\link{dtwclust}}.
#' Note that for multivariate series, this means that it \strong{must} be a list of matrices, even if the list
#' has only one element.
#'
#' @exportMethod predict
#'
setMethod("predict", "dtwclust",
          function(object, newdata = NULL, ...) {
               if (is.null(newdata)) {
                    if (object@type != "fuzzy")
                         ret <- object@cluster
                    else
                         ret <- object@fcluster

               } else {
                    newdata <- consistency_check(newdata, "tsmat")
                    nm <- names(newdata)
                    consistency_check(newdata, "vltslist")
                    newdata <- object@family@preproc(newdata)
                    distmat <- object@family@dist(newdata, object@centers, ...)
                    ret <- object@family@cluster(distmat = distmat, m = object@control@fuzziness)

                    if (object@type != "fuzzy")
                         names(ret) <- nm
                    else
                         dimnames(ret) <- list(nm, paste0("cluster_", 1:ncol(ret)))
               }

               ret
          })

# ========================================================================================================
# Plot
# ========================================================================================================

#' @details
#' The plot method uses the \code{ggplot2} plotting system (see \code{\link[ggplot2]{ggplot}}).
#'
#' The default depends on whether a hierarchical method was used or not. In those cases, the dendrogram is
#' plotted by default; you can pass any extra parameters to \code{\link[stats]{plot.hclust}} via \code{...}.
#'
#' Otherwise, the function plots the time series of each cluster along with the obtained centroid.
#' The default values for cluster centroids are: \code{linetype = "dashed"}, \code{size = 1.5},
#' \code{colour = "black"}, \code{alpha = 0.5}. You can change this by means of \code{...}.
#'
#' You can choose what to plot with the \code{type} parameter. Possible options are:
#'
#' \itemize{
#'   \item \code{"dendrogram"}: Only available for hierarchical clustering.
#'   \item \code{"series"}: Plot the time series divided into clusters without including centroids.
#'   \item \code{"centroids"}: Plot the obtained centroids only.
#'   \item \code{"sc"}: Plot both series and centroids
#' }
#'
#' The flag \code{save.data} should be set to \code{TRUE} when running \code{\link{dtwclust}} to be able to
#' use this. Optionally, you can manually provide the data in the \code{data} parameter.
#'
#' If created, the function returns the \code{gg} object invisibly, in case you want to modify it to your
#' liking. You might want to look at \code{\link[ggplot2]{ggplot_build}} if that's the case.
#'
#' @param y Ignored.
#' @param ... Further arguments to pass to \code{\link[ggplot2]{geom_line}} for the plotting of the
#' \emph{cluster centers}, or to \code{\link[stats]{plot.hclust}}. See details.
#' @param clus A numeric vector indicating which clusters to plot.
#' @param labs.arg Arguments to change the title and/or axis labels. See \code{\link[ggplot2]{labs}} for more
#' information
#' @param data The data in the same format as it was provided to \code{\link{dtwclust}}.
#' @param time Optional values for the time axis. If series have different lengths, provide the time values of
#' the longest series.
#' @param plot Logical flag. You can set this to \code{FALSE} in case you want to save the ggplot object without
#' printing anything to screen
#' @param type What to plot. \code{NULL} means default. See details.
#' @param show.centroids Deprecated.
#'
#' @return The plot method returns a \code{gg} object (or \code{NULL} for dendrogram plot) invisibly.
#'
#' @rdname dtwclust-methods
#' @aliases plot,dtwclust,missing-method
#'
#' @exportMethod plot
#'
setMethod("plot", signature(x = "dtwclust", y = "missing"),
          function(x, y, ...,
                   clus = seq_len(x@k), labs.arg = NULL,
                   data = NULL, time = NULL,
                   plot = TRUE, type = NULL,
                   show.centroids = TRUE) {

               if (!missing(show.centroids))
                    warning("The 'show.centroids' argument has been deprecated. Use 'type' instead.")

               ## set default type if none was provided
               if (!is.null(type))
                    type <- match.arg(type, c("dendrogram", "series", "centroids", "sc"))
               else if (x@type == "hierarchical")
                    type <- "dendrogram"
               else
                    type <- "sc"

               ## plot dendrogram?
               if(x@type == "hierarchical" && type == "dendrogram") {
                    x <- S3Part(x, strictS3 = TRUE)
                    if (plot) plot(x, ...)
                    return(invisible(NULL))

               } else if (x@type != "hierarchical" && type == "dendrogram") {
                    ## dendrogram for non-hierarchical is not possible
                    stop("Dendrogram plot only applies to hierarchical clustering.")
               }

               ## Obtain data, the priority is: provided data > included data list
               if (!is.null(data)) {
                    data <- consistency_check(data, "tsmat")

                    Lengths <- lengths(data)
                    L <- max(Lengths)
                    trail <- L - Lengths

                    df <- mapply(data, trail, SIMPLIFY = FALSE,
                                 FUN = function(series, trail) {
                                      c(series, rep(NA, trail))
                                 })

               } else if (length(x@datalist) > 0){
                    Lengths <- lengths(x@datalist)
                    L <- max(Lengths)
                    trail <- L - Lengths

                    df <- mapply(x@datalist, trail, SIMPLIFY = FALSE,
                                 FUN = function(series, trail) {
                                      c(series, rep(NA, trail))
                                 })
               } else {
                    stop("Provided object has no data. Please re-run the algorithm with save.data = TRUE
                         or provide the data manually.")
               }

               ## Obtain centers (which can be matrix or lists of series)
               Lengths <- lengths(x@centers)
               trail <- L - Lengths

               cen <- mapply(x@centers, trail, SIMPLIFY = FALSE,
                             FUN = function(series, trail) {
                                  c(series, rep(NA, trail))
                             })

               cen <- as.data.frame(cen)
               colnames(cen) <- NULL

               ## Check if data was z-normalized
               if (x@preproc == "zscore")
                    titleStr <- "Clusters' members (z-normalized)"
               else
                    titleStr <- "Clusters' members"

               ## transform data

               df <- as.data.frame(df)

               if (is.null(time)) {
                    t <- seq_len(L)

               } else {
                    if (length(time) != L)
                         stop("Length mismatch between values and time stamps")

                    t <- time
               }

               df <- data.frame(t = t, df)
               dfm <- reshape2::melt(df, id.vars = "t")

               cl <- rep(x@cluster, each = L)
               color <- lapply(tabulate(x@cluster), function(i) {
                    rep(1L:i, each = L)
               })
               color <- factor(unlist(color))

               dfm <- data.frame(dfm, cl = cl, color = color)

               ## transform centers

               cen <- data.frame(t = t, cen)
               cenm <- reshape2::melt(cen, id.vars = "t")
               cl <- rep(1L:x@k, each = L)
               cenm <- data.frame(cenm, cl = cl)

               ## create gg object
               gg <- ggplot(data.frame(t = integer(),
                                       variable = factor(),
                                       value = numeric(),
                                       cl = factor(),
                                       color = factor()),
                            aes_string(x = "t",
                                       y = "value",
                                       group = "variable"))

               if (type %in% c("sc", "centroids")) {
                    if(length(list(...)) == 0L)
                         gg <- gg + geom_line(data = cenm[cenm$cl %in% clus, ],
                                              linetype = "dashed",
                                              size = 1.5,
                                              colour = "black",
                                              alpha = 0.5)
                    else
                         gg <- gg + geom_line(data = cenm[cenm$cl %in% clus, ], ...)
               }

               if (type %in% c("sc", "series")) {
                    gg <- gg + geom_line(data = dfm[dfm$cl %in% clus, ], aes(colour = color))
               }

               gg <- gg +
                    facet_wrap(~cl, scales = "free_y") +
                    guides(colour = FALSE) +
                    theme_bw()

               if (!is.null(labs.arg))
                    gg <- gg + labs(labs.arg)
               else
                    gg <- gg + labs(title = titleStr)

               ## If NAs were introduced by me (for length consistency), I want them to be removed
               ## automatically, and I don't want a warning
               if (plot) suppressWarnings(plot(gg))

               invisible(gg)
          })

# ========================================================================================================
# Cluster validity indices
# ========================================================================================================

#' @title
#' Cluster validity indices
#'
#' @description
#' Compute different cluster validity indices (CVIs) of a given cluster partition, using the
#' clustering distance measure if applicable.
#'
#' @details
#'
#' Clustering is commonly considered to be an unsupervised procedure, so evaluating its performance
#' can be rather subjective. However, a great amount of effort has been invested in trying to standardize
#' cluster evaluation metrics by using cluster validity indices (CVIs).
#'
#' CVIs can be classified as internal, external or relative depending on how they are computed.
#' Focusing on the first two, the crucial difference is that internal CVIs only consider the partitioned
#' data and try to define a measure of cluster purity, whereas external CVIs compare the obtained partition
#' to the correct one. Thus, external CVIs can only be used if the ground truth is known. Each index defines
#' their range of values and whether they are to be minimized or maximized. In many cases, these CVIs can be
#' used to evaluate the result of a clustering algorithm regardless of how the clustering works internally,
#' or how the partition came to be.
#'
#' Knowing which CVI will work best cannot be determined a priori, so they should be tested for each
#' specific application. Usually, many CVIs are utilized and compared to each other, maybe using a majority
#' vote to decide on a final result. Furthermore, it should be noted that many CVIs perform additional
#' distance calculations when being computed, which can be very considerable if using DTW.
#'
#' Note that, even though a fuzzy partition can be changed into a crisp one, making it compatible with many
#' of the existing CVIs, there are also fuzzy CVIs tailored specifically to fuzzy clustering, and these may
#' be more suitable in those situations, but have not been implemented here yet.
#'
#' External CVIs (the first 4 are calculated via \code{\link[flexclust]{comPart}}):
#'
#' \itemize{
#'   \item \code{"RI"}: Rand Index (to be maximized).
#'   \item \code{"ARI"}: Adjusted Rand Index (to be maximized).
#'   \item \code{"J"}: Jaccard Index (to be maximized).
#'   \item \code{"FM"}: Fowlkes-Mallows (to be maximized).
#'   \item \code{"VI"}: Variation of Information (to be minimized).
#' }
#'
#' Internal CVIs:
#'
#' \itemize{
#'   \item \code{"Sil"}: Silhouette index (to be maximized). A normalized summation-type index based on
#'   intracluster distance and the distance to nearest neighbors. This index essentially calculates (or
#'   re-uses, if already available) the whole distance matrix using the clustering distance measure, so if
#'   you were trying to avoid that in the first place, this may not be the best CVI for your application.
#' }
#'
#' @name cvi
#' @rdname cvi
#'
#' @param a An object returned by the \code{\link{dtwclust}} function, or a vector that can be coerced to
#' integers which indicate the cluster memeberships.
#' @param b If needed, a vector that can be coerced to integers which indicate the cluster memeberships.
#' @param type Character vector indicating which indices are to be computed. See details.
#' @param ... Arguments to pass to and from other methods.
#' @param log.base Base of the logarithm to be used in the calculation of VI.
#'
#' @return The chosen CVIs
#'
#' @references
#'
#' Arbelaitz, O., Gurrutxaga, I., Muguerza, J., Perez, J. M., & Perona, I. (2013). An extensive
#' comparative study of cluster validity indices. Pattern Recognition, 46(1), 243-256.
#'
#' @exportMethod cvi
#'
setGeneric("cvi", def = function(a, b = NULL, type = "valid", ..., log.base = 2) {
     ## Only external CVIs here
     if (is.null(b))
          stop("A second set of cluster membership indices is required in 'b' for this/these CVI(s).")

     a <- as.integer(a)
     b <- as.integer(b)

     type <- match.arg(type, several.ok = TRUE,
                       c("RI", "ARI", "J", "FM", "VI",
                         "valid", "external"))

     if (any(type %in% c("valid", "external")))
          type <- c("RI", "ARI", "J", "FM", "VI")

     which_flexclust <- type %in% c("RI", "ARI", "J", "FM")

     if (any(which_flexclust))
          CVIs <- flexclust::comPart(x = a, y = b, type = type[which_flexclust])
     else
          CVIs <- numeric()

     if (any(type == "VI")) {
          ## Variation of information
          ## taken from https://github.com/cran/mcclust/blob/master/R/vi.dist.R

          ## entropy
          ent <- function(cl) {
               n <- length(cl)
               p <- table(cl) / n
               -sum(p * log(p, base = log.base))
          }

          ## mutual information
          mi <- function(cl1, cl2) {
               p12 <- table(cl1, cl2) / length(cl1)
               p1p2 <- outer(table(cl1) / length(cl1), table(cl2) / length(cl2))
               sum(p12[p12 > 0] * log(p12[p12 > 0] / p1p2[p12 > 0], base = log.base))
          }

          VI <- ent(a) + ent(b) - 2 * mi(a, b)
          CVIs <- c(CVIs, VI = VI)
     }

     CVIs
})

#' @rdname cvi
#' @aliases cvi,dtwclust-method
#'
setMethod("cvi", signature(a = "dtwclust"),
          function(a, b = NULL, type = "valid", ...) {
               if (a@type == "fuzzy")
                    stop("Only CVIs for crisp partitions are currently implemented.")

               type <- match.arg(type, several.ok = TRUE,
                                 c("RI", "ARI", "J", "FM", "VI",
                                   "Sil", "SF", "CH", "DB", "DB*", "DB**", "D",
                                   "valid", "internal", "external"))

               dots <- list(...)

               internal <- c("Sil", "SF", "CH", "DB", "DB*", "DB**", "D")
               external <- c("RI", "ARI", "J", "FM", "VI")

               if (any(type == "valid")) {
                    type <- if(is.null(b)) internal else c(internal, external)

               } else if (any(type == "internal")) {
                    type <- internal

               } else if (any(type == "external")) {
                    type <- external
               }

               which_internal <- type %in% internal
               which_external <- type %in% external

               if (any(which_external))
                    CVIs <- cvi(a@cluster, b = b, type = type[which_external], ...)
               else
                    CVIs <- numeric()

               if (any(which_internal)) {
                    if (is.null(a@distmat)) {
                         distmat <- do.call(a@family@dist,
                                            args = c(list(x = a@datalist, centers = NULL),
                                                     a@dots))
                    } else {
                         distmat <- a@distmat
                    }

                    CVIs <- c(CVIs, sapply(type[which_internal], function(CVI) {
                         switch(EXPR = CVI,
                                ## Silhouette index
                                Sil = {
                                     c_k <- as.numeric(table(a@cluster)[a@cluster])

                                     ab <- lapply(unique(a@cluster), function(k) {
                                          idx <- a@cluster == k

                                          this_a <- rowSums(distmat[idx, idx, drop = FALSE]) / c_k[idx]

                                          this_b <- apply(distmat[idx, !idx, drop = FALSE], 1L, function(row) {
                                               ret <- row / c_k[!idx]
                                               ret <- min(tapply(ret, a@cluster[!idx], sum))
                                               ret
                                          })

                                          data.frame(a = this_a, b = this_b)
                                     })

                                     ab <- do.call(rbind, ab)

                                     Sil <- sum((ab$b - ab$a) / apply(ab, 1L, max)) / length(a@datalist)

                                     Sil
                                },

                                ## Default for now
                                -1)
                    }))
               }

               CVIs
          })

# ========================================================================================================
# Rand Index from flexclust package
# ========================================================================================================

#' @title
#' Compare partitions
#'
#' @description
#' Compute the (adjusted) Rand, Jaccard and Fowlkes-Mallows index for agreement of two partitions.
#' This generic is included in the \code{flexclust} package.
#'
#' @seealso
#'
#' \code{\link[flexclust]{randIndex}}
#'
#' @name randIndex
#' @rdname randIndex
#'
#' @param x,y,correct,original See \code{\link[flexclust]{randIndex}}.
#'
#' @exportMethod randIndex
#'
NULL

#' @rdname randIndex
#' @aliases randIndex,dtwclust,ANY-method
#'
setMethod("randIndex", signature(x="dtwclust", y="ANY"),
          function(x, y, correct = TRUE, original = !correct) {
               randIndex(x@cluster, y, correct = correct, original = original)
          })

#' @rdname randIndex
#' @aliases randIndex,ANY,dtwclust-method
#'
setMethod("randIndex", signature(x="ANY", y="dtwclust"),
          function(x, y, correct = TRUE, original = !correct) {
               randIndex(x, y@cluster, correct = correct, original = original)
          })

#' @rdname randIndex
#' @aliases randIndex,dtwclust,dtwclust-method
#'
setMethod("randIndex", signature(x="dtwclust", y="dtwclust"),
          function(x, y, correct = TRUE, original = !correct) {
               randIndex(x@cluster, y@cluster, correct = correct, original = original)
          })

# ========================================================================================================
# Cluster Similarity from flexclust
# ========================================================================================================

#' @title
#' Cluster Similarity Matrix
#'
#' @description
#' Returns a matrix of cluster similarities. Currently two methods for computing
#' similarities of clusters are implemented.
#' This generic is included in the \code{flexclust} package.
#'
#' @seealso
#'
#' \code{\link[flexclust]{clusterSim}}
#'
#' @name clusterSim
#' @rdname clusterSim
#'
#' @param object,data,method,symmetric,... See \code{\link[flexclust]{clusterSim}}.
#'
#' @exportMethod clusterSim
#'
#' @rdname clusterSim
#' @aliases clusterSim,dtwclust-method
#'
setMethod("clusterSim", "dtwclust",
          function (object, data = NULL,
                    method = c("shadow", "centers"),
                    symmetric = FALSE, ...)
          {
               method <- match.arg(method)

               if (object@k == 1)
                    return(matrix(1))

               if (method == "shadow") {
                    if (is.null(data))
                         data <- object@datalist
                    else
                         data <- consistency_check(data, "tsmat")

                    distmat <- object@family@dist(data, object@centers)

                    r <- t(matrix(apply(distmat, 1, rank, ties.method = "first"),
                                  nrow = ncol(distmat)))
                    cluster <- lapply(1:2, function(k) {
                         apply(r, 1, function(x) which(x == k))
                    })

                    K <- max(cluster[[1]])
                    z <- matrix(0, ncol = K, nrow = K)
                    for (k in 1:K) {
                         ok1 <- cluster[[1]] == k
                         if (any(ok1)) {
                              for (n in 1:K) {
                                   if (k != n) {
                                        ok2 <- ok1 & cluster[[2]] == n
                                        if (any(ok2)) {
                                             z[k, n] <- 2 * sum(distmat[ok2, k] / (distmat[ok2, k] + distmat[ok2, n]))
                                        }
                                   }
                              }
                              z[k, ] <- z[k, ]/sum(ok1)
                         }
                    }

                    diag(z) <- 1

                    if (symmetric)
                         z <- (z + t(z))/2

               } else {
                    z <- object@family@dist(object@centers, object@centers)
                    z <- 1 - z/max(z)
               }

               z
          })

# ========================================================================================================
# Validity and coercion methods for control
# ========================================================================================================

setValidity("dtwclustControl",
            function(object) {

                 if (!is.null(object@window.size) && object@window.size < 1)
                      return("Window size must be positive if provided")

                 object@norm <- match.arg(object@norm, c("L1", "L2"))

                 if (object@dba.iter < 0L)
                      return("DBA iterations must be positive")

                 if (object@iter.max < 0L)
                      return("Maximum iterations must be positive")

                 if (object@nrep < 1L)
                      return("Number of repetitions must be at least one")

                 if (object@fuzziness <= 1)
                      return("Fuzziness exponent should be greater than one")

                 if (object@delta < 0)
                      return("Delta should be positive")

                 TRUE
            })

setAs("list", "dtwclustControl",
      function(from, to) {
           ctrl <- new(to)

           num <- c("delta", "fuzziness")

           for (arg in names(from)) {
                val <- from[[arg]]

                if (is.numeric(val) && !(arg %in% num))
                     val <- as.integer(val)

                slot(ctrl, arg) <- val
           }

           validObject(ctrl)

           ctrl
      })

setAs("NULL", "dtwclustControl",
      function(from, to) {
           new(to)
      })

# ========================================================================================================
# Coercion methods for cross/pair-dist
# ========================================================================================================

as.matrix.crossdist <- function(x, ...) {
     class(x) <- NULL
     as.matrix(x, ...)
}

as.matrix.pairdist <- function(x, ...) {
     class(x) <- NULL
     as.matrix(x, ...)
}

as.data.frame.crossdist <- function(x, ...) {
     as.data.frame(as.matrix(x, ...), ...)
}

as.data.frame.pairdist <- function(x, ...) {
     as.data.frame(as.matrix(x, ...), ...)
}
