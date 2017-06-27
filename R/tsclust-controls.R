#' Control parameters for clusterings with [tsclust()]
#'
#' Control parameters for fine-grained control.
#'
#' @name tsclust-controls
#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param pam.precompute Logical flag. Precompute the whole distance matrix once and reuse it on
#'   each iteration if using PAM centroids. Otherwise calculate distances at every iteration. See
#'   details.
#' @param iter.max Integer. Maximum number of allowed iterations for partitional/fuzzy clustering.
#' @param nrep Integer. How many times to repeat clustering with different starting points.
#' @param symmetric Logical flag. Is the distance function symmetric? In other words, is `dist(x,y)`
#'   == `dist(y,x)`? If `TRUE`, only half the distance matrix needs to be computed. Automatically
#'   detected and overridden for the distances included in \pkg{dtwclust}.
#' @param packages Character vector with the names of any packages required for custom `proxy`
#'   functions. Since the distance entries are re-registered in each parallel worker if needed, this
#'   is probably useless, but just in case.
#' @param distmat If available, the cross-distance matrix can be provided here. Only relevant for
#'   partitional with PAM centroids or hierarchical procedures.
#' @param pam.sparse Attempt to use a sparse matrix for PAM centroids. See details.
#'
#' @details
#'
#' The functions essentially return their function arguments in a classed list, although some checks
#' are performed.
#'
#' @section Partitional:
#'
#'   As of version 4.0.0, using `pam.precompute = FALSE` is optimized by defining a sparse matrix
#'   (see [Matrix::sparseMatrix()]) and updating it every iteration (except for `"dtw_lb"`
#'   distance). For smaller datasets, precomputing the whole distance matrix is still probably
#'   faster.
#'
#'   Explicitly setting both `pam.precompute` and `pam.sparse` to `FALSE` will reproduce previous
#'   behavior, and might be faster if the distance matrix is very big but the distance function is
#'   very quick.
#'
#'   Parallel computations for PAM centroids have the following considerations:
#'
#'   - If `pam.precompute` is `TRUE`, both distance matrix calculations and repetitions are done in
#'     parallel, regardless of `pam.sparse`.
#'   - If `pam.precompute` is `FALSE` and `pam.sparse` is `TRUE`, repetitions are done sequentially,
#'     so that the distance calculations can be done in parallel and the sparse matrix updated
#'     iteratively.
#'   - If both `pam.precompute` and `pam.sparse` are `FALSE`, repetitions are done in parallel, and
#'     each repetition performs distance calculations sequentially, but the sparse matrix cannot be
#'     updated iteratively.
#'
partitional_control <- function(pam.precompute = TRUE,
                                iter.max = 100L,
                                nrep = 1L,
                                symmetric = FALSE,
                                packages = character(0L),
                                distmat = NULL,
                                pam.sparse = !pam.precompute)
{
    if (any(iter.max <= 0L)) stop("Maximum iterations must be positive")
    if (any(nrep < 1L)) stop("Number of repetitions must be at least one")

    structure(
        list(pam.precompute = as.logical(pam.precompute),
             pam.sparse = as.logical(pam.sparse),
             iter.max = as.integer(iter.max),
             nrep = as.integer(nrep)[1L],
             symmetric = as.logical(symmetric)[1L],
             packages = unique(c("dtwclust", as.character(packages))),
             distmat = distmat),
        "class" = c(control_classes[["partitional"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param method Character vector with one or more linkage methods to use in hierarchical procedures
#'   (see [stats::hclust()]), the character `"all"` to use all of the available ones, or a function
#'   that performs hierarchical clustering based on distance matrices (e.g. [cluster::diana()]). See
#'   details.
#'
#' @section Hierarchical:
#'
#'   There are some limitations when using a custom hierarchical function in `method`: it will
#'   receive the lower triangular of the distance matrix as first argument (see [stats::as.dist()])
#'   and the result should support the [stats::as.hclust()] generic. This functionality was added
#'   with the \pkg{cluster} package in mind, since its functions follow this convention, but other
#'   functions could be used if they are adapted to work similarly.
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

        if ("all" %in% method)
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
        "class" = c(control_classes[["hierarchical"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param fuzziness Numeric. Exponent used for fuzzy clustering. Commonly termed `m` in the
#'   literature.
#' @param delta Numeric. Convergence criterion for fuzzy clustering.
#'
fuzzy_control <- function(fuzziness = 2,
                          iter.max = 100L,
                          delta = 1e-3,
                          packages = character(0L),
                          symmetric = FALSE)
{
    if (any(fuzziness <= 1)) stop("Fuzziness exponent should be greater than one")
    if (any(iter.max <= 0L)) stop("Maximum iterations must be positive")
    if (any(delta < 0)) stop("Delta should be positive")

    structure(
        list(fuzziness = fuzziness,
             iter.max = as.integer(iter.max),
             delta = delta,
             symmetric = as.logical(symmetric)[1L],
             packages = unique(c("dtwclust", as.character(packages)))),
        "class" = c(control_classes[["fuzzy"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param dc The cutoff distance for the TADPole algorithm.
#' @param window.size The window.size specifically for the TADPole algorithm.
#' @param lb The lower bound to use with TADPole. Either `"lbk"` or `"lbi"`.
#'
tadpole_control <- function(dc,
                            window.size,
                            lb = "lbk")
{
    if (any(dc <= 0)) stop("Cutoff distance 'dc' must be positive")
    window.size <- check_consistency(window.size, "window")
    lb <- match.arg(lb, c("lbk", "lbi"), several.ok = TRUE)

    structure(
        list(dc = dc,
             window.size = window.size,
             lb = lb),
        "class" = c(control_classes[["tadpole"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param preproc A list of arguments for a preprocessing function to be used in [tsclust()].
#' @param dist A list of arguments for a distance function to be used in [tsclust()].
#' @param cent A list of arguments for a centroid function to be used in [tsclust()].
#'
#' @section TADPole:
#'
#'   When using TADPole, the `dist` argument list includes the `window.size` and specifies `norm =
#'   "L2"`.
#'
tsclust_args <- function(preproc = list(), dist = list(), cent = list())
{
    structure(
        list(preproc = preproc,
             dist = dist,
             cent = cent),
        "class" = c(control_classes[["args"]])
    )
}

