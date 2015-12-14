#' @title Methods for \code{dtwclust}
#'
#' @description
#' Methods associated with \code{\link{dtwclust-class}} objects.
#'
#' @name dtwclust-methods
#' @rdname dtwclust-methods
#'
#' @seealso \code{\link{dtwclust-class}}, \code{\link{dtwclust}}, \code{\link[ggplot2]{ggplot}}
#'
NULL

#' @details
#' The plot method plots the time series of each cluster along with the obtained centroid.
#' It uses \code{ggplot2} plotting system (\code{\link[ggplot2]{ggplot}}).
#'
#' Note that if \code{type} was \code{"hierarchical"} when \code{\link{dtwclust}} was called, the dendrogram
#' will be plotted instead, and no object returned.
#'
#' The flag \code{save.data} should be set to \code{TRUE} when running \code{\link{dtwclust}} to be able to
#' use this. Optionally, you can manually provide the data in the \code{data} parameter.
#'
#' The function returns the \code{gg} object invisibly, in case you want to modify it to your liking. You
#' might want to look at \code{\link[ggplot2]{ggplot_build}} if that's the case.
#'
#' @param x An object of class \code{\link{dtwclust-class}} as returned by \code{\link{dtwclust}}.
#' @param y Ignored.
#' @param clus A numeric vector indicating which clusters to plot.
#' @param labs.arg Arguments to change the title and/or axis labels. See \code{\link[ggplot2]{labs}} for more
#' information
#' @param data The data in the same format as it was provided to \code{\link{dtwclust}}.
#' @param time Optional values for the time axis. If series have different lengths, provide the time values of
#' the longest series.
#' @param plot Logical flag. You can set this to \code{FALSE} in case you want to save the ggplot object without
#' printing anything to screen
#' @param ... Further arguments to pass to \code{\link[ggplot2]{geom_line}} for the plotting of the
#' \emph{cluster centers}. Default values are: \code{linetype = "dashed"}, \code{size = 1.5},
#' \code{colour = "black"}, \code{alpha = 0.5}.
#'
#' @return The plot method returns a \code{gg} object (or \code{NULL} for hierarchical methods) invisibly.
#'
#' @rdname dtwclust-methods
#' @aliases plot,dtwclust,missing-method
#'
#' @exportMethod plot
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
setMethod("plot", signature(x="dtwclust", y="missing"),
          function(x, y, ..., clus = seq_len(x@k),
                   labs.arg = NULL, data = NULL, time = NULL, plot = TRUE) {

               ## If .Data part is not empty, object has 'hclust' part
               if(length(x@.Data) > 0) {
                    x <- S3Part(x, strictS3 = TRUE)
                    plot(x)
                    return(invisible(NULL))
               }

               ## Obtain data, the priority is: provided data > included data matrix > included data list
               if (!is.null(data)) {
                    if (is.matrix(data))
                         data <- consistency_check(data, "tsmat")

                    lengths <- sapply(data, length)
                    L <- max(lengths)
                    trail <- L - lengths

                    df <- mapply(data, trail, SIMPLIFY = FALSE,
                                 FUN = function(series, trail) {
                                      c(series, rep(NA, trail))
                                 })

               } else if (length(x@datalist) > 0){
                    lengths <- sapply(x@datalist, length)
                    L <- max(lengths)
                    trail <- L - lengths

                    df <- mapply(x@datalist, trail, SIMPLIFY = FALSE,
                                 FUN = function(series, trail) {
                                      c(series, rep(NA, trail))
                                 })
               } else {
                    stop("Provided object has no data. Please re-run the algorithm with save.data = TRUE
                         or provide the data manually.")
               }

               ## Obtain centers (which can be matrix or lists of series)
               lengths <- sapply(x@centers, length)
               trail <- L - lengths

               cen <- mapply(x@centers, trail, SIMPLIFY = FALSE,
                             FUN = function(series, trail) {
                                  c(series, rep(NA, trail))
                             })

               cen <- as.data.frame(cen)

               ## Check if data was z-normalized
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

               df <- cbind(t, df)
               dfm <- reshape2::melt(df, id.vars = "t")

               cl <- rep(x@cluster, each = L)
               color <- sapply(tabulate(x@cluster), function(i) {
                    rep(1:i, each = L)
               })
               color <- factor(unlist(color))

               dfm <- cbind(dfm, cl, color)

               ## tarnsform centers

               cen <- cbind(t, cen)
               cenm <- reshape2::melt(cen, id.vars = "t")
               cl <- rep(1:x@k, each = L)
               cenm <- cbind(cenm, cl)

               ## create gg object

               if (length(list(...)) == 0) {
                    g <- ggplot(dfm[dfm$cl %in% clus, ], aes_string(x="t", y="value", group="variable")) +
                         geom_line(data = cenm[cenm$cl %in% clus, ], linetype = "dashed", size = 1.5,
                                   colour = "black", alpha = 0.5) +
                         geom_line(aes(colour = color)) +
                         facet_wrap(~cl, scales = "free_y") +
                         guides(colour=FALSE) +
                         theme_bw()

               } else {
                    g <- ggplot(dfm[dfm$cl %in% clus, ], aes_string(x="t", y="value", group="variable")) +
                         geom_line(data = cenm[cenm$cl %in% clus, ], ...) +
                         geom_line(aes(colour = color)) +
                         facet_wrap(~cl, scales = "free_y") +
                         guides(colour=FALSE) +
                         theme_bw()

               }

               if (!is.null(labs.arg))
                    g <- g + labs(labs.arg)
               else
                    g <- g + labs(title = titleStr)

               ## If NAs were introduced by me (for length consistency), I want them to be removed
               ## automatically, and I don't want a warning
               if (plot)
                    suppressWarnings(print(g))

               invisible(g)
          })


#' @details
#' Show method displays basic information of results.
#'
#' @rdname dtwclust-methods
#' @aliases show,dtwclust-method
#'
#' @param object An object of class \code{dtwclust}.
#'
setMethod("show", "dtwclust",
          function(object) {
               cat(object@type, "clustering with", object@k, "clusters\n")
               cat("Using", object@distance, "distance\n")

               if (object@type == "partitional")
                    cat("Using", object@centroid, "centroids\n")
               if (object@preproc != "none")
                    cat("Using", object@preproc, "preprocessing\n")

               cat("\nTime required for analysis:\n")
               print(object@proctime)

               cat("\nCluster sizes with average intra-cluster distance:\n\n")
               print(object@clusinfo)
          })


#' Compare partitions
#'
#' Compute the (adjusted) Rand, Jaccard and Fowlkes-Mallows index for agreement of two partitions.
#'
#' This generic is included in the \code{flexclust} package. It will no longer be included with \code{dtwclust}
#' in the next release.
#'
#' @seealso
#'
#' \code{\link[flexclust]{randIndex}}
#'
#' @name randIndex
#' @rdname randIndex
#' @aliases randIndex,dtwclust,ANY-method
#'
#' @exportMethod randIndex
#'
#' @param x,y,correct,original See \code{\link[flexclust]{randIndex}}.
#'
#' @return A vector of indices.
#'
setMethod("randIndex", signature(x="dtwclust", y="ANY"),
          function(x, y, correct = TRUE, original = !correct) {
               warning("Package 'dtwclust': Support for this generic method will be dropped in the next release.")

               callNextMethod()
          })

#' @rdname randIndex
#' @aliases randIndex,ANY,dtwclust-method
#'
setMethod("randIndex", signature(x="ANY", y="dtwclust"),
          function(x, y, correct = TRUE, original = !correct) {
               warning("Package 'dtwclust': Support for this generic method will be dropped in the next release.")

               callNextMethod()
          })

#' @rdname randIndex
#' @aliases randIndex,dtwclust,dtwclust-method
#'
setMethod("randIndex", signature(x="dtwclust", y="dtwclust"),
          function(x, y, correct = TRUE, original = !correct) {
               warning("Package 'dtwclust': Support for this generic method will be dropped in the next release.")

               callNextMethod()
          })
