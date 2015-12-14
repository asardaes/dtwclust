#' Class definition for \code{dtwclust}
#'
#' Formal S4 class to know how to handle data for plotting.
#'
#' The class will no longer inherit from \code{\link[flexclust]{kccasimple-class}} in the next release. However,
#' it now contains \code{hclust} as superclass.
#'
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

setClass("dtwclust", contains = c("kccasimple", "hclust"),
         slots = c(k="integer",
                   cluster="integer",
                   iter="integer",
                   converged="logical",
                   clusinfo="data.frame",

                   centers="list",
                   cldist="matrix",

                   type = "character",
                   distance = "character",
                   centroid = "character",
                   preproc = "character",
                   datalist = "list",
                   proctime = "proc_time"))
