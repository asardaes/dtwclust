#' Control parameters for clusterings with \code{\link{tsclust}}
#'
#' Control parameters for fine-grained control.
#'
#' @name tsclust-controls
#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param pam.precompute Logical flag. Precompute the whole distance matrix once and reuse it on
#'   each iteration if using PAM centroids. Otherwise calculate distances at every iteration.
#' @param iter.max Integer. Maximum number of allowed iterations for partitional/fuzzy clustering.
#' @param nrep Integer. How many times to repeat clustering with different starting points.
#' @param symmetric Logical flag. Is the distance function symmetric? In other words, is
#'   \code{dist(x,y)} == \code{dist(y,x)}? If \code{TRUE}, only half the distance matrix needs to be
#'   computed. Overridden if the function detects an invalid user-provided value.
#' @param packages Character vector with the names of any packages required for custom \code{proxy}
#'   functions. Since the distance entries are re-registered in each parallel worker if needed, this
#'   is probably useless, but just in case.
#' @param distmat If available, the cross-distance matrix can be provided here. Only relevant for
#'   partitional with PAM centroids or hierarchical procedures.
#'
#' @details
#'
#' The functions essentially return their function arguments in a classed list, although some checks
#' are performed.
#'
partitional_control <- function(pam.precompute = TRUE,
                                iter.max = 100L,
                                nrep = 1L,
                                symmetric = FALSE,
                                packages = character(0L),
                                distmat = NULL)
{
    if (iter.max <= 0L) stop("Maximum iterations must be positive")

    if (nrep < 1L) stop("Number of repetitions must be at least one")

    structure(
        list(pam.precompute = as.logical(pam.precompute),
             iter.max = as.integer(iter.max),
             nrep = as.integer(nrep),
             symmetric = as.logical(symmetric),
             packages = unique(c("dtwclust", as.character(packages))),
             distmat = distmat),
        "class" = c("PtCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param method Character vector with one or more linkage methods to use in hierarchical procedures
#'   (see \code{\link[stats]{hclust}}), the character \code{"all"} to use all of the available ones,
#'   or a function that performs hierarchical clustering based on distance matrices (e.g.
#'   \code{\link[cluster]{diana}}).
#'
hierarchical_control <- function(method = "average",
                                 symmetric = FALSE,
                                 packages = character(0L),
                                 distmat = NULL)
{
    if (is.character(method)) {
        method <- match.arg(method,
                            c("ward.D", "ward.D2", "single", "complete",
                              "average", "mcquitty", "median", "centroid",
                              "all"),
                            several.ok = TRUE)

        if (any(method == "all"))
            method <- c("ward.D", "ward.D2", "single", "complete",
                        "average", "mcquitty", "median", "centroid")

    } else if (!is.function(method))
        stop("Argument 'method' must be either a supported character or a function.")
    else
        attr(method, "name") <- as.character(substitute(method))[1L]

    structure(
        list(method = method,
             symmetric = as.logical(symmetric),
             packages = unique(c("dtwclust", as.character(packages))),
             distmat = distmat),
        "class" = c("HcCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param fuzziness Numeric. Exponent used for fuzzy clustering. Commonly termed \code{m} in the
#'   literature.
#' @param delta Numeric. Convergence criterion for fuzzy clustering.
#'
fuzzy_control <- function(fuzziness = 2,
                          iter.max = 100L,
                          delta = 1e-3,
                          packages = character(0L))
{
    if (fuzziness <= 1) stop("Fuzziness exponent should be greater than one")

    if (iter.max <= 0L) stop("Maximum iterations must be positive")

    if (delta < 0) stop("Delta should be positive")

    structure(
        list(fuzziness = fuzziness,
             iter.max = as.integer(iter.max),
             delta = delta,
             packages = unique(c("dtwclust", as.character(packages)))),
        "class" = c("FzCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param dc The cutoff distance for the TADPole algorithm.
#' @param window.size The window.size specifically for the TADPole algorithm.
#' @param lb The lower bound to use with TADPole. Either \code{"lbk"} or \code{"lbi"}.
#'
tadpole_control <- function(dc,
                            window.size,
                            lb = "lbk")
{
    if (dc <= 0) stop("Cutoff distance 'dc' must be positive")

    window.size <- check_consistency(window.size, "window")

    lb <- match.arg(lb, c("lbk", "lbi"))

    structure(
        list(dc = dc,
             window.size = window.size,
             lb = lb),
        "class" = c("TpCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param preproc A list of arguments for a preprocessing function to be used in
#'   \code{\link{tsclust}}.
#' @param dist A list of arguments for a distance function to be used in \code{\link{tsclust}}.
#' @param cent A list of arguments for a centroid function to be used in \code{\link{tsclust}}.
#'
#' @details
#'
#' When using TADPole, the \code{dist} argument list includes the \code{window.size} and specifies
#' \code{norm = "L2"}.
#'
tsclust_args <- function(preproc = list(), dist = list(), cent = list())
{
    structure(
        list(preproc = preproc,
             dist = dist,
             cent = cent),
        "class" = c("TscArgs")
    )
}

