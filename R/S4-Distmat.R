# ==================================================================================================
# Distmat RC and methods to transparently handle PAM centroids
# ==================================================================================================

#' Distance matrix
#'
#' Reference class that is used internally for cross-distance matrices.
#'
#' @importFrom methods is
#' @importFrom methods setRefClass
#'
#' @field distmat A distance matrix.
#' @field series Time series list.
#' @field distfun The distance function to calculate the distance matrix.
#' @field dist_args Arguments for the distance function.
#' @field id_cent Indices of the centroids (if any).
#'
#' @keywords internal
#'
Distmat <- methods::setRefClass("Distmat",
                       fields = list(
                           distmat = "ANY",
                           series = "list",
                           distfun = "function",
                           dist_args = "list",
                           id_cent = "integer"
                       ),
                       methods = list(
                           initialize = function(..., distmat, series, distance, control, error.check = TRUE) {
                               "Initialization based on needed parameters"

                               if (missing(distmat)) {
                                   if (tolower(distance) == "dtw_lb") distance <- "dtw_basic"
                                   if (error.check) {
                                       check_consistency(series, "vltslist")
                                       check_consistency(distance,
                                                         "dist",
                                                         trace = FALSE,
                                                         diff_lengths = different_lengths(series),
                                                         silent = FALSE)
                                       if (!methods::is(control, "PtCtrl"))
                                           stop("Invalid control provided.") # nocov
                                   }
                                   # need another dist closure, otherwise it would be recursive
                                   control$distmat <- NULL
                                   initFields(...,
                                              series = series,
                                              distfun = ddist2(distance, control))
                               }
                               else
                                   initFields(..., distmat = distmat)
                               # return
                               invisible(NULL)
                           }
                       )
)

#' Generics for `Distmat`
#'
#' Generics with methods for [Distmat-class].
#'
#' @name Distmat-generics
#' @rdname Distmat-generics
#' @keywords internal
#' @importFrom methods setMethod
#'
NULL

#' @rdname Distmat-generics
#' @aliases [,Distmat,ANY,ANY,ANY
#'
#' @param x A [Distmat-class] object.
#' @param i Row indices.
#' @param j Column indices.
#' @param ... Ignored.
#' @param drop Logical to drop dimensions after subsetting.
#'
#' @details
#'
#' Accessing matrix elements with `[]` first calculates the values if necessary.
#'
setMethod(`[`, "Distmat", function(x, i, j, ..., drop = TRUE) {
    if (inherits(x$distmat, "uninitializedField")) {
        if (inherits(x$distfun, "uninitializedField"))
            stop("Invalid internal Distmat instance.") # nocov

        centroids <- if (identical(i,j)) NULL else x$series[j]
        dm <- quoted_call(x$distfun, x = x$series[i], centroids = centroids, dots = x$dist_args)
    }
    else {
        dm <- x$distmat[i, j, drop = drop]
        if (identical(dim(dm), dim(x$distmat))) attributes(dm) <- attributes(x$distmat)
    }
    # return
    dm
})

dim.Distmat <- function(x) { dim(x$distmat) } # nocov
