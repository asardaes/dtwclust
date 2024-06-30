lower_triangular_index <- function(i, j, n, diagonal) {
    stopifnot(i > 0L, i <= n, j > 0L, j <= i)

    i <- i - 1L
    j <- j - 1L

    adjustment <- if (diagonal) {
        0L
    }
    else {
        Reduce(x = 0L:j, init = 0L, f = function(a, b) { a + b + 1L })
    }

    i + j * n - adjustment + 1L
}

#' Distance matrix's lower triangular
#'
#' Reference class that is used internally for PAM centroids when `pam.precompute = TRUE` and
#' `pam.sparse = FALSE`. It contains [Distmat-class].
#'
#' @include S4-Distmat.R
#' @importFrom methods setRefClass
#'
#' @field distmat The lower triangular.
#'
#' @keywords internal
#'
DistmatLowerTriangular <- methods::setRefClass(
    "DistmatLowerTriangular",
    contains = "Distmat",
    fields = list(
        distmat = "ANY"
    ),
    methods = list(
        initialize = function(..., distmat) {
            "Initialization based on needed parameters"

            if (missing(distmat)) {
                stop("distmat must be provided for this class.")
            }
            else if (!inherits(distmat, "dist")) {
                stop("distmat must be a 'dist' object.")
            }

            callSuper(..., distmat = distmat)
            # return
            invisible(NULL)
        }
    )
)

#' Generics for `DistmatLowerTriangular`
#'
#' Generics with methods for [DistmatLowerTriangular-class].
#'
#' @name DistmatLowerTriangular-generics
#' @rdname DistmatLowerTriangular-generics
#' @keywords internal
#' @importFrom methods setMethod
#'
NULL

#' @rdname DistmatLowerTriangular-generics
#' @aliases show,DistmatLowerTriangular
#' @importFrom methods show
#'
#' @param object A [DistmatLowerTriangular-class] object.
#'
setMethod("show", "DistmatLowerTriangular", function(object) { methods::show(object$distmat) }) # nocov

#' @rdname DistmatLowerTriangular-generics
#' @aliases [,DistmatLowerTriangular,ANY,ANY,ANY
#'
#' @param x A [DistmatLowerTriangular-class] object.
#' @param i Row indices.
#' @param j Column indices.
#' @param ... Ignored.
#'
setMethod(`[`, "DistmatLowerTriangular", function(x, i, j, ...) {
    if (missing(j)) {
        stopifnot(inherits(i, "matrix"), ncol(i) == 2L)
        j <- i[, 2L]
        i <- i[, 1L]
        drop <- TRUE
    }
    else {
        out_dim <- c(length(i), length(j))
        out_dimnames <- list(i, j)
        combinations <- expand.grid(i = i, j = j)
        i <- combinations$i
        j <- combinations$j
        drop <- FALSE

    }

    n <- attr(x$distmat, "Size")
    diagonal <- isTRUE(attr(x$distmat, "Diag"))
    entries <- mapply(i, j, FUN = function(i, j) {
        if (!diagonal && i == j) {
            0
        }
        else if (j > i) {
            x$distmat[lower_triangular_index(j, i, n, diagonal)]
        }
        else {
            x$distmat[lower_triangular_index(i, j, n, diagonal)]
        }
    })

    if (drop) {
        entries
    }
    else {
        dim(entries) <- out_dim
        dimnames(entries) <- out_dimnames
        entries
    }
})

#' @exportS3Method base::dim
dim.DistmatLowerTriangular <- function(x) { rep(attr(x$distmat, "Size"), 2L) } # nocov

methods::setOldClass("dist")
methods::setOldClass("crossdist")

setAs("matrix", "Distmat", function(from) { Distmat$new(distmat = from) })
setAs("crossdist", "Distmat", function(from) { Distmat$new(distmat = from) })
setAs("dist", "Distmat", function(from) { DistmatLowerTriangular$new(distmat = from) })

#' @exportS3Method base::as.matrix
as.matrix.distdiag <- function(x) {
    n <- attr(x, "Size")
    m <- matrix(0, n, n)
    m[lower.tri(m, diag = TRUE)] <- x

    lbls <- attr(x, "Labels")
    if (!is.null(lbls)) {
        names(m) <- list(lbls, lbls)
    }

    m
}
