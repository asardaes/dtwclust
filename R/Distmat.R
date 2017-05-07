# ==================================================================================================
# Distmat RC and methods to transparently handle PAM's pam.precompute = FALSE non-sparse case
# ==================================================================================================

#' Distance matrix
#'
#' Reference class that is used internally for PAM centroids when both `pam.precompute` and
#' `pam.sparse` are `FALSE`.
#'
#' @field distmat A distance matrix.
#' @field series Time series list.
#' @field distfun The distance function to calculate the distance matrix.
#' @field dist_args Arguments for the distance function.
#' @field id_cent Indices of the centroids (if any).
#'
Distmat <- setRefClass("Distmat",
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
                                   if (error.check) {
                                       check_consistency(series, "vltslist")
                                       check_consistency(distance,
                                                         "dist",
                                                         trace = FALSE,
                                                         Lengths = different_lengths(series),
                                                         silent = FALSE)

                                       if (class(control) != "PtCtrl")
                                           stop("Invalid control provided.")
                                   }

                                   ## need another dist closure, otherwise it would be recursive
                                   control$distmat <- NULL
                                   initFields(...,
                                              series = series,
                                              distfun = ddist2(distance, control))

                               } else
                                   initFields(..., distmat = distmat)

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
    dm_exists <- !inherits(x$distmat, "uninitializedField")

    if (!dm_exists) {
        dm <- do.call(x$distfun,
                      enlist(x = x$series[i],
                             centroids = x$series[j],
                             dots = x$dist_args))
    } else {
        dm <- x$distmat[i, j, drop = drop]
        if (identical(dim(dm), dim(x$distmat))) attributes(dm) <- attributes(x$distmat)
    }

    dm
})

dim.Distmat <- function(x) { dim(x$distmat) }
