#' Class definition for \code{dtwclust}
#'
#' Formal S4 class to know how to handle data for plotting.
#'
#' The class inherits from \code{\link[flexclust]{kccasimple-class}}, so most related slots and
#' methods are also supported.
#'
#' @slot type A string indicating one of the supported clustering types of \code{\link{dtwclust}}.
#' @slot distance A string indicating the distance used with \code{\link{dtwclust}}.
#' @slot centroid A string indicating the centroid used with \code{\link{dtwclust}}.
#' @slot preproc A string indicating the preprocessing used with \code{\link{dtwclust}}.
#' @slot datalist The provided data in the form of a list, where each element is a time series.
#' @slot proctime Time during function execution, as measured by \code{\link[base]{proc.time}}.
#'
#' @name dtwclust-class
#' @rdname dtwclust-class
#' @import methods
#' @exportClass dtwclust
#'
NULL

setClass("proc_time4", contains = "numeric", slots = c(names = "character"))
setOldClass("proc_time", S4Class = "proc_time4")
removeClass("proc_time4")

setClass("dtwclust", contains = c("kccasimple"),
         slots = c(type = "character",
                   distance = "character",
                   centroid = "character",
                   preproc = "character",
                   datalist = "list",
                   proctime = "proc_time"))
