#' Class definition for \code{dtwclust}
#'
#' Formal S4 class to know how to handle data for plotting.
#'
#' It contains the following specific slots:
#'
#' \itemize{
#'   \item \code{type}: A string indicating one of the supported clustering types of \code{\link{dtwclust}}.
#'   \item \code{distance}: A string indicating the distance used with \code{\link{dtwclust}}.
#'   \item \code{centroid}: A string indicating the centroid used with \code{\link{dtwclust}}.
#'   \item \code{preproc}: A string indicating the preprocessing used with \code{\link{dtwclust}}.
#'   \item \code{datalist}: The provided data in the form of a list, where each element is a time series.
#' }
#'
#' Additionally, the class inherits from \code{\link[flexclust]{kccasimple-class}}, so most related slots and
#' methods are also supported.
#'
#' @name dtwclust-class
#' @rdname dtwclust-class
#' @import methods
#' @exportClass dtwclust
#'

setClass("dtwclust", contains = c("kccasimple"),
         slots = c(type = "character",
                   distance = "character",
                   centroid = "character",
                   preproc = "character",
                   datalist = "list"))
