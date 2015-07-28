#' Class definition for \code{dtwclust}
#'
#' Formal S4 class to know how to handle data for plotting.
#'
#' It contains the following specific slots:
#'
#' \itemize{
#'   \item \code{type}: A string indicating one of the supported clustering types of \code{\link{dtwclust}}.
#'   \item \code{distance}: A string indicating one of the supported distances of \code{\link{dtwclust}}.
#'   \item \code{centroid}: A string indicating one of the supported centroids of \code{\link{dtwclust}}.
#' }
#'
#' Additionally, the class inherits from \code{\link[flexclust]{kccasimple-class}}, so all related slots and
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
                   centroid = "character"))
