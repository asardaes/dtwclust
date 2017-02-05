#' Control parameters for clusterings with \code{\link{tsclust}}
#'
#' Control parameters
#'
#' @name tsclust-controls
#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @details
#'
#' Partitional
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
             packages = c("dtwclust", as.character(packages)),
             distmat = distmat),
        "class" = c("list", "PtCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @details
#'
#' Hierarchical
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
             packages = c("dtwclust", as.character(packages)),
             distmat = distmat),
        "class" = c("list", "HcCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @details
#'
#' Fuzzy
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
             packages = c("dtwclust", as.character(packages))),
        "class" = c("list", "FzCtrl")
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @details
#'
#' TADPole
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
        "class" = c("list", "TpCtrl")
    )
}
