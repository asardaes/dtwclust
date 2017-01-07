#' Time series warping envelops
#'
#' This function computes the envelops for DTW lower bound calculations with a Sakoe-Chiba band for
#' a given univariate time series using the streaming algorithm proposed by Lemire (2009).
#'
#' @export
#'
#' @param x A univariate time series.
#' @param window.size Window size for envelop calculation. See details.
#' @param error.check Check data inconsistencies?
#'
#' @details
#'
#' The windowing constraint uses a centered window. The calculations expect a value in
#' \code{window.size} that represents the distance between the point considered and one of the edges
#' of the window. Therefore, if, for example, \code{window.size = 10}, the warping for an
#' observation \eqn{x_i} considers the points between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting
#' in \code{10(2) + 1 = 21} observations falling within the window.
#'
#' @return A list with two elements (lower and upper envelops): \code{lower} and \code{upper}.
#'
#' @note
#'
#' This envelop is calculated assuming a Sakoe-Chiba constraint for DTW.
#'
#' @references
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .''
#' \emph{Pattern Recognition}, \strong{42}(9), pp. 2169 - 2180. ISSN 0031-3203,
#' \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030},
#' \url{http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
#'
compute_envelop <- function(x, window.size, error.check = TRUE) {
    if (error.check) {
        if (is_multivariate(list(x)))
            stop("The envelop can conly be computed for univariate series.")

        check_consistency(x, "ts")
    }

    window.size <- check_consistency(window.size, "window")
    window.size <- window.size * 2L + 1L

    ## NOTE: window.size in this function is window.size*2 + 1, thus the 2L below
    if (window.size > (2L * length(x)))
        stop("Window cannot be greater or equal than the series' length.")

    .Call("envelop", x, window.size, PACKAGE = "dtwclust")
}
