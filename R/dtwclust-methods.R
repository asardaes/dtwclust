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
#' The plot method, by default, plots the time series of each cluster along with the obtained centroid.
#' It uses \code{ggplot2} plotting system (see \code{\link[ggplot2]{ggplot}}). The default values for
#' cluster centers are: \code{linetype = "dashed"}, \code{size = 1.5}, \code{colour = "black"},
#' \code{alpha = 0.5}. You can change this by means of \code{...}.
#'
#' The flag \code{save.data} should be set to \code{TRUE} when running \code{\link{dtwclust}} to be able to
#' use this. Optionally, you can manually provide the data in the \code{data} parameter.
#'
#' The function returns the \code{gg} object invisibly, in case you want to modify it to your liking. You
#' might want to look at \code{\link[ggplot2]{ggplot_build}} if that's the case.
#'
#' If a hierarchical procedure was used, then you can specify \code{type} \code{=} \code{"dendrogram"} to
#' plot the corresponding dendrogram (the default in this case), and pass any extra parameters via \code{...}.
#' Use \code{type} \code{=} \code{"series"} to plot the time series clusters using the original call's \code{k}.
#'
#' @param y Ignored.
#' @param ... Further arguments to pass to \code{\link[ggplot2]{geom_line}} for the plotting of the
#' \emph{cluster centers}, or to \code{\link[stats]{plot.hclust}}. See details.
#' @param clus A numeric vector indicating which clusters to plot.
#' @param labs.arg Arguments to change the title and/or axis labels. See \code{\link[ggplot2]{labs}} for more
#' information
#' @param show.centroids Logical flag. Should cluster centroids be included in the plots?
#' @param data The data in the same format as it was provided to \code{\link{dtwclust}}.
#' @param time Optional values for the time axis. If series have different lengths, provide the time values of
#' the longest series.
#' @param plot Logical flag. You can set this to \code{FALSE} in case you want to save the ggplot object without
#' printing anything to screen
#' @param type What to plot. Only relevant for hierarchical procedures. See details.
#'
#' @return The plot method returns a \code{gg} object (or \code{NULL} for hierarchical methods) invisibly.
#'
#' @rdname dtwclust-methods
#' @aliases plot,dtwclust,missing-method
#'
#' @exportMethod plot
#'
setMethod("plot", signature(x="dtwclust", y="missing"),
          function(x, y, ..., clus = seq_len(x@k),
                   labs.arg = NULL, show.centroids = TRUE,
                   data = NULL, time = NULL,
                   plot = TRUE, type = "dendrogram") {

               type <- match.arg(type, c("dendrogram", "series"))

               ## plot dendrogram?
               if(x@type == "hierarchical" && type == "dendrogram") {
                    x <- S3Part(x, strictS3 = TRUE)
                    if (plot) plot(x, ...)
                    return(invisible(NULL))
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
                    titleStr <- "Clusters and their members, including cluster center (all z-normalized)"
               else
                    titleStr <- "Clusters and their members, including cluster center"

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
               gg <- ggplot(dfm[dfm$cl %in% clus, ], aes_string(x = "t", y = "value", group = "variable"))

               if (show.centroids) {
                    if(length(list(...)) == 0L)
                         gg <- gg + geom_line(data = cenm[cenm$cl %in% clus, ],
                                              linetype = "dashed",
                                              size = 1.5,
                                              colour = "black",
                                              alpha = 0.5)
                    else
                         gg <- gg + geom_line(data = cenm[cenm$cl %in% clus, ], ...)
               }

               gg <- gg +
                    geom_line(aes(colour = color)) +
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
# info from modeltools (perhaps for internal support of some flexclust functions)
# ========================================================================================================

# setMethod("info", "dtwclust",
#           function(object, which, ...) {
#                ret <- switch(EXPR = which,
#                              help = c("distsum", "size", "av_dist"),
#                              distsum = sum(object@cldist[, 1]),
#                              size = {
#                                   var <- object@clusinfo$size
#                                   names(var) <- rownames(object@clusinfo)
#
#                                   var
#                              },
#                              av_dist = {
#                                   var <- object@clusinfo$av_dist
#                                   names(var) <- rownames(object@clusinfo)
#
#                                   var
#                              },
#
#                              stop("Requested info not available. Use which = 'help'."))
#
#                ret
#           })

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
