# ==================================================================================================
# Sparse distmat RC and methods to transparently handle PAM's pam.precompute = FALSE case
# ==================================================================================================

#' Sparse distance matrix
#'
#' Reference class that is used internally for PAM centroids when `pam.precompute = FALSE`. It
#' allows for mutable state. It contains [Distmat-class].
#'
#' @include Distmat.R
#' @include pkg.R
#'
#' @field distmat The sparse matrix.
#' @field symmetric Logical indicating if the matrix is symmetric.
#' @field distmat_indices External pointer (C++ class) with the indices of existing values within
#'   the matrix, and the method to update them.
#'
SparseDistmat <- setRefClass(
    "SparseDistmat",
    contains = "Distmat",
    fields = list(
        distmat = "sparseMatrix",
        symmetric = "logical",
        distmat_indices = "externalptr"
    ),
    methods = list(
        initialize = function(..., control) {
            "Initialization based on needed parameters"

            callSuper(..., control = control)
            symmetric <<- control$symmetric

            distmat <<- Matrix::sparseMatrix(i = 1L:length(series),
                                             j = 1L:length(series),
                                             x = 0,
                                             symmetric = control$symmetric)

            if (isTRUE(control$symmetric) && distmat@uplo != "L") distmat <<- t(distmat)

            ## initialize C++ class
            distmat_indices <<- .Call(C_SparseDistmatIndices__new,
                                      nrow(distmat),
                                      PACKAGE = "dtwclust")
            sapply(1L:length(series), function(id) {
                .Call(C_SparseDistmatIndices__getNewIndices,
                      distmat_indices, id, id, symmetric,
                      PACKAGE = "dtwclust")
            })

            invisible(NULL)
        }
    )
)

#' Generics for `SparseDistmat`
#'
#' Generics with methods for [SparseDistmat-class].
#'
#' @name SparseDistmat-generics
#' @rdname SparseDistmat-generics
#'
NULL

#' @rdname SparseDistmat-generics
#' @aliases show,SparseDistmat
#'
#' @param object A [SparseDistmat-class] object.
#'
setMethod("show", "SparseDistmat", function(object) { show(object$distmat) }) # nocov

#' @rdname SparseDistmat-generics
#' @aliases [,SparseDistmat,ANY,ANY,ANY
#'
#' @param x A [SparseDistmat-class] object.
#' @param i Row indices.
#' @param j Column indices.
#' @param ... Ignored.
#' @param drop Logical to drop dimensions after subsetting.
#'
#' @details
#'
#' Accessing matrix elements with `[]` first calculates the values if necessary.
#'
setMethod(`[`, "SparseDistmat", function(x, i, j, ..., drop = TRUE) {
    id_new <- .Call(C_SparseDistmatIndices__getNewIndices,
                    x$distmat_indices, i, j, x$symmetric,
                    PACKAGE = "dtwclust")

    ## update distmat if necessary
    if (nrow(id_new) > 0L) {
        x$distmat[id_new] <- as.numeric(do.call(x$distfun,
                                                enlist(x = x$series[id_new[, 1L]],
                                                       centroids = x$series[id_new[, 2L]],
                                                       pairwise = TRUE,
                                                       dots = x$dist_args)))

        if (x$symmetric) x$distmat <- Matrix::forceSymmetric(x$distmat, "L")
    }

    x$distmat[i, j, drop = drop]
})

dim.SparseDistmat <- function(x) { dim(x$distmat) } # nocov
