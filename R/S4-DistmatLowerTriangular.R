lower_triangular_index <- function(i, j, n, diagonal) {
    stopifnot(i > 0L, i <= n, j > 0L, j <= i)
    if (!diagonal) stopifnot(i != j)

    diagonal <- as.integer(diagonal)
    adjustment <- Reduce(x = 1L:j, init = 0L, f = function(a, b) {
        a + b - diagonal
    })

    i + (j - 1L) * n - adjustment
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
#' @details
#'
#' If you wish to, you can use this class to access `dist` elements with `[` as if it were a normal
#' matrix. You can use [methods::new] passing the `dist` object in a `distmat` argument.
#'
#' @examples
#'
#' dm <- new("DistmatLowerTriangular",
#'           distmat = proxy::dist(CharTraj[1:5], method = "gak", sigma = 5.5, window.size = 10L))
#'
#' dm[2:3, 4:5]
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

lti <- function(i) {
    if (is.logical(i)) which(i) else i
}

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
        stopifnot(inherits(i, "matrix"))
        if (is.logical(i)) {
            stopifnot(attr(x$distmat, "Size") == dim(i))
            j <- rep(1:ncol(i), apply(i, 2L, sum))
            i <- as.integer(apply(i, 2L, which))
        }
        else {
            stopifnot(ncol(i) == 2L)
            j <- lti(i[, 2L])
            i <- lti(i[, 1L])
        }
        drop <- TRUE
    }
    else {
        i <- lti(i)
        if (any(i < 0L)) {
            i <- seq_len(attr(x$distmat, "Size"))[i]
        }

        j <- lti(j)
        if (any(j < 0L)) {
            j <- seq_len(attr(x$distmat, "Size"))[j]
        }

        out_dim <- c(length(i), length(j))
        out_dimnames <- list(i, j)
        combinations <- expand.grid(i = i, j = j)
        i <- combinations$i
        j <- combinations$j
        drop <- FALSE
    }

    n <- attr(x$distmat, "Size")
    diagonal <- attr(x$distmat, "method") == "SDTW" & isTRUE(attr(x$distmat, "Diag"))
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
