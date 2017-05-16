# ==================================================================================================
# Sparse distmat RC and methods to transparently handle PAM's pam.precompute = FALSE case
# ==================================================================================================

#' Sparse distance matrix
#'
#' Reference class that is used internally for PAM centroids when `pam.precompute = FALSE`. It
#' allows for mutable state. It contains [Distmat-class].
#'
#' @include Distmat.R
#'
#' @field distmat The sparse matrix.
#' @field symmetric Logical indicating if the matrix is symmetric.
#' @field existing_ids Matrix with the indices of existing values within the matrix.
#'
SparseDistmat <- setRefClass("SparseDistmat",
                             contains = "Distmat",
                             fields = list(
                                 distmat = "sparseMatrix",
                                 symmetric = "logical",
                                 existing_ids = "matrix"
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

                                     if (isTRUE(control$symmetric) && distmat@uplo != "L")
                                         distmat <<- t(distmat)

                                     existing_ids <<- base::as.matrix(
                                         Matrix::summary(distmat)[c("i", "j")]
                                     )

                                     invisible(NULL)
                                 }
                             ))

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
setMethod("show", "SparseDistmat", function(object) { show(object$distmat) })

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
    ## number of rows of existing indices
    rows <- nrow(x$existing_ids)
    ## indices of needed vals
    id_new <- base::as.matrix(expand.grid(i = i, j = j))

    ## modify indices to lower triangular if necessary (symmetric only)
    if (x$symmetric)
        id_new <- t(apply(id_new, 1L, function(this_row) {
            if (this_row[2L] > this_row[1L]) this_row[2L:1L] else this_row
        }))

    ## extract only indices of needed values that don't exist yet
    id_new <- rbind(x$existing_ids, id_new)
    id_duplicated <- duplicated(id_new)
    rows <- (rows + 1L):nrow(id_new)
    id_new <- id_new[rows, , drop = FALSE][!id_duplicated[rows], , drop = FALSE]

    ## update distmat if necessary
    if (nrow(id_new) > 0L) {
        x$distmat[id_new] <- as.numeric(do.call(x$distfun,
                                                enlist(x = x$series[id_new[, 1L]],
                                                       centroids = x$series[id_new[, 2L]],
                                                       pairwise = TRUE,
                                                       dots = x$dist_args)))

        if (x$symmetric) x$distmat <- Matrix::forceSymmetric(x$distmat, "L")
        x$existing_ids <- rbind(x$existing_ids, id_new)
    }

    x$distmat[i, j, drop = drop]
})

dim.SparseDistmat <- function(x) { dim(x$distmat) }
