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
#' @param labs.arg Arguments to change the title and/or axis labels. See \code{\link[ggplot2]{labs}} for more
#' information
#' @param data The data in the same format as it was provided to \code{\link{dtwclust}}.
#' @param time Optional values for the time axis. If series have different lengths, provide the time values of
#' the longest series.
#' @param plot Boolean flag. You can set this to FALSE in case you want to save the ggplot object without
#' printing anything to screen
#' @param ... Further arguments to pass to \code{\link[ggplot2]{geom_line}} for the plotting of the
#' \emph{cluster centers}. Default values are: \code{linetype = "dashed"}, \code{size = 1.5},
#' \code{colour = "black"}, \code{alpha = 0.5}.
#'
#' @exportMethod plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom modeltools empty
#'
NULL


#' @rdname plot-methods
#' @aliases plot,dtwclust,missing-method
#'
setMethod("plot", signature(x="dtwclust", y="missing"),
          function(x, y, clus=seq_len(x@k),
                   labs.arg = NULL, data=NULL, time = NULL, plot = TRUE, ...) {

               ## Obtain data, the priority is: provided data > included data matrix > included data list

               if (!is.null(data)) {
                    if (is.matrix(data))
                         df <- t(data)

                    else if (is.list(data)) {
                         lengths <- sapply(data, length)
                         L <- max(lengths)
                         trail <- L - lengths

                         df <- mapply(data, trail, SIMPLIFY = FALSE,
                                      FUN = function(series, trail) {
                                           c(series, rep(NA, trail))
                                      })
                    }

               } else if (!modeltools::empty(x@data)) {
                    df <- t(x@data@get("designMatrix"))
                    L <- nrow(df)

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

               df <- as.data.frame(df)

               ## Obtain centers (which can be matrix or lists of series)

               if (is.matrix(x@centers)) {
                    cen <- as.data.frame(t(x@centers))

               } else if (is.list(x@centers)) {
                    lengths <- sapply(x@centers, length)
                    trail <- L - lengths

                    cen <- mapply(x@centers, trail, SIMPLIFY = FALSE,
                                  FUN = function(series, trail) {
                                       c(series, rep(NA, trail))
                                  })

                    cen <- as.data.frame(cen)
               }

               ## Check if data was z-normalized

               if (x@preproc == "zscore") {
                    df <- apply(df, 2, zscore)
                    titleStr <- "Clusters and their members, including cluster center (all z-normalized)"
               } else {
                    titleStr <- "Clusters and their members, including cluster center"
               }

               ## transform data

               #n <- nrow(df)
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
