#' Centroid for partition around medoids
#'
#' Extract the medoid time series based on a distance measure.
#'
#' @export
#' @importFrom Matrix rowSums
#'
#' @param series The time series in one of the formats accepted by [tslist()].
#' @param distance A character indicating which distance to use. Only needed if `distmat` is `NULL`.
#'   The distance must be registered in [proxy::pr_DB()].
#' @param ids Integer vector indicating which of the `series` should be considered.
#' @param distmat Optionally, a pre-computed cross-distance matrix of *all* `series`.
#' @param ... Any extra parameters for the `distance` function that may be used.
#' @template error-check
#'
#' @details
#'
#' The medoid's index is determined by taking the \eqn{arg min} of the `distmat`'s row-sums
#' (considering only the rows in `ids`). The distance matrix is calculated if needed.
#'
#' @return
#'
#' The medoid time series.
#'
#' @examples
#'
#' # Computes the distance matrix for all series
#' pam_cent(CharTraj, "dtw_basic", ids = 6L:10L, window.size = 15L) # series_id = 7L
#'
#' # Computes the distance matrix for the chosen subset only
#' pam_cent(CharTraj[6L:10L], "dtw_basic", window.size = 15L) # series_id = 2L
#'
pam_cent <- function(series, distance, ids = seq_along(series), distmat = NULL, ...,
                     error.check = TRUE)
{
    series <- tslist(series, simplify = error.check)
    dots <- list(...)
    i_cl <- dots$.i_cl_
    dots$.i_cl_ <- NULL

    if (!inherits(distmat, "Distmat")) {
        if (missing(distance))
            distance <- attr(distmat, "method")

        args <- list(
            distmat = distmat,
            series = series,
            dist_args = dots,
            distance = distance,
            control = partitional_control(),
            error.check = error.check
        )

        if (is.null(distmat)) {
            if (is.null(distance))
                stop("If 'distmat' is missing, 'distance' must be provided.")

            args$distmat <- NULL
        }

        # S4-Distmat.R
        distmat <- do.call(Distmat$new, args, TRUE)
    }

    d <- distmat[ids, ids, drop = FALSE]
    d <- rowSums(d) # d can be normal matrix or from Matrix package, so no namespace here
    id_cent <- ids[which.min(d)]
    cent <- series[[id_cent]]

    if (is.null(i_cl))
        attr(cent, "series_id") <- id_cent
    else
        distmat$id_cent[i_cl] <- id_cent

    cent
}
