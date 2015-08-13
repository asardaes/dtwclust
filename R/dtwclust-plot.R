#' Plot the result of \code{dtwclust}
#'
#' Plots the time series of each cluster along with the obtained centroid. It uses \code{ggplot2} plotting
#' system.
#'
#' The flag \code{save.data} must be set to \code{TRUE} when running \code{\link{dtwclust}} to be able to
#' use this.
#'
#' Optionally, you can manually provide the clustering result as well as the data in \code{data}.
#'
#' The function returns the \code{gg} object invisibly, in case you want to modify it to your liking. You
#' might want to look at \code{\link[ggplot2]{ggplot_build}} if that's the case.
#'
#' @name plot-dtwclust
#' @docType methods
#' @rdname plot-methods
#'
#' @seealso \code{\link{dtwclust-class}}, \code{\link{dtwclust}}, \code{\link[ggplot2]{ggplot}}
#'
#' @param x An object of class \code{\link{dtwclust-class}} as returned by \code{\link{dtwclust}}.
#' @param y Ignored.
#' @param clus Which clusters to plot.
#' @param data The data in the same format as it was provided to \code{\link{dtwclust}}.
#' @param labs.arg Arguments to change the title and/or axis labels. See \code{\link[ggplot2]{labs}} for more
#' information
#' @param ... Further arguments to pass to \code{\link[ggplot2]{geom_line}} for the plotting of the
#' \emph{cluster centers}. Default values are: \code{linetype = "dashed"}, \code{size = 1.5},
#' \code{colour = "black"}, \code{alpha = 0.5}.
#'
#' @exportMethod plot
#' @import ggplot2
#' @importFrom reshape2 melt
#'
NULL


#' @rdname plot-methods
#' @aliases plot,dtwclust,missing-method
#'
setMethod("plot", signature(x="dtwclust", y="missing"),
          function(x, y, clus=seq_len(x@k), data=NULL, labs.arg = NULL, ...) {

               if (!is.null(data))
                    df <- t(data)
               else if (!is.null(x@data))
                    df <- t(x@data@get("designMatrix"))
               else
                    stop("Provided object hast no data. Please re-run the algorithm with save.data = TRUE
                         or provide the data manually.")

               if (x@centroid == "shape") {
                    df <- apply(df, 2, zscore)
                    titleStr <- "Clusters and their members, including cluster center (all z-normalized)"
               } else {
                    titleStr <- "Clusters and their members, including cluster center"
               }

               df <- as.data.frame(df)

               n <- nrow(df)
               t <- seq_len(n)
               df <- cbind(t, df)
               dfm <- melt(df, id.vars = "t")

               cl <- rep(x@cluster, each = n)
               color <- sapply(tabulate(x@cluster), function(i) {
                    rep(1:i, each = n)
               })
               color <- factor(unlist(color))

               dfm <- cbind(dfm, cl, color)

               cen <- as.data.frame(t(x@centers))
               cen <- cbind(t, cen)
               cenm <- melt(cen, id.vars = "t")
               cl <- rep(1:x@k, each = n)
               cenm <- cbind(cenm, cl)

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

               print(g)

               invisible(g)
          })
