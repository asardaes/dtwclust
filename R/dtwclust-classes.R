#' Class definition for \code{dtwclustFamily}
#'
#' Formal S4 class with the family of functions used for partitional and hierarchical procedures in
#' \code{\link{dtwclust}}.
#'
#' The custom implementations also handle parallelization.
#'
#' @slot dist The function to calculate the distance matrices.
#' @slot allcent The function to calculate centroids at each iteration.
#' @slot cluster The function used to assign a series to a cluster.
#' @slot preproc The function used to preprocess the data (important for \code{\link[stats]{predict}}).
#'
#' @name dtwclustFamily-class
#' @rdname dtwclustFamily-class
#' @aliases dtwclustFamily
#'
#' @import methods
#' @exportClass dtwclustFamily
#'
setClass("dtwclustFamily",

         slots = c(dist = "function",
                   allcent = "function",
                   cluster = "function",
                   preproc = "function"),

         prototype = prototype(preproc = function(x, ...) x,
                               cluster = function(x, centers, distmat = NULL) {
                                    if (is.null(distmat))
                                         distmat <- dist(x, centers)

                                    max.col(-distmat, "first")
                               })
)

#' Class definition for \code{dtwclust}
#'
#' Formal S4 class to know how to handle data for plotting.
#'
#' The class no longer inherits from \code{\link[flexclust]{kccasimple-class}}. However,
#' it now contains \code{hclust} as superclass, and most slots were ported. Namely, \code{data} slot wasn't.
#'
#' @slot call The function call.
#' @slot family An object of class \code{\link{dtwclustFamily}}.
#' @slot k Integer indicating the number of desired clusters.
#' @slot cluster Integer vector indicating which cluster a series belongs to.
#' @slot iter The number of iterations used.
#' @slot converged A logical indicating whether the function converged.
#' @slot clusinfo A data frame with two columns: \code{size} indicates the number of series each cluster has, and
#' \code{av_dist} indicates the average distance between series for each cluster.
#' @slot centers A list with the centroid time series.
#' @slot cldist A column vector with the distance between each series in the data and its corresponding centroid.
#' @slot type A string indicating one of the supported clustering types of \code{\link{dtwclust}}.
#' @slot distance A string indicating the distance used with \code{\link{dtwclust}}.
#' @slot centroid A string indicating the centroid used with \code{\link{dtwclust}}.
#' @slot preproc A string indicating the preprocessing used with \code{\link{dtwclust}}.
#' @slot datalist The provided data in the form of a list, where each element is a time series.
#' @slot proctime Time during function execution, as measured by \code{\link[base]{proc.time}}.
#'
#' @name dtwclust-class
#' @rdname dtwclust-class
#'
#' @import methods
#' @exportClass dtwclust
#'
NULL

setClass("proc_time4", contains = "numeric", slots = c(names = "character"))
setOldClass("proc_time", S4Class = "proc_time4")
removeClass("proc_time4")

setClass("hclust4", contains = "list", slots = c(names = "character"))
setOldClass("hclust", S4Class = "hclust4")
removeClass("hclust4")

setClass("dtwclust", contains = c("hclust"),
         slots = c(call = "call",
                   family = "dtwclustFamily",

                   k = "integer",
                   cluster = "integer",
                   iter = "integer",
                   converged = "logical",
                   clusinfo = "data.frame",

                   centers = "list",
                   cldist = "matrix",

                   type = "character",
                   distance = "character",
                   centroid = "character",
                   preproc = "character",
                   datalist = "list",
                   proctime = "proc_time"))
