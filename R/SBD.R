#' Shape-based distance
#'
#' Distance based on coefficient-normalized cross-correlation as proposed by Paparrizos and Gravano
#' (2015) for the k-Shape clustering algorithm.
#'
#' @export
#'
#' @param x,y Univariate time series.
#' @param znorm Logical. Should each series be z-normalized before calculating the distance?
#' @template error-check
#'
#' @details
#'
#' This distance works best if the series are *z-normalized*. If not, at least they should have
#' appropriate amplitudes, since the values of the signals **do** affect the outcome.
#'
#' If `x` and `y` do **not** have the same length, it would be best if the longer sequence is
#' provided in `y`, because it will be shifted to match `x`. After matching, the series may have to
#' be truncated or extended and padded with zeros if needed.
#'
#' The output values lie between 0 and 2, with 0 indicating perfect similarity.
#'
#' @return A list with:
#'
#'   - `dist`: The shape-based distance between `x` and `y`.
#'   - `yshift`: A shifted version of `y` so that it optimally matches `x` (based on [NCCc()]).
#'
#' @template proxy
#' @template symmetric
#' @section Proxy version:
#'
#'   In some situations, e.g. for relatively small distance matrices, the overhead introduced by the
#'   logic that computes only half the distance matrix can be bigger than just calculating the whole
#'   matrix.
#'
#' @note
#'
#' If you wish to calculate the distance between several time series, it would be better to use the
#' version registered with the `proxy` package, since it includes some small optimizations. See the
#' examples.
#'
#' This distance is calculated with help of the Fast Fourier Transform, so it can be sensitive to
#' numerical precision. Thus, this function (and the functions that depend on it) might return
#' different values in 32 bit installations compared to 64 bit ones.
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In *Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data*, series
#' SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' @seealso
#'
#' [NCCc()], [shape_extraction()]
#'
#' @examples
#'
#' # load data
#' data(uciCT)
#'
#' # distance between series of different lengths
#' sbd <- SBD(CharTraj[[1]], CharTraj[[100]], znorm = TRUE)$dist
#'
#' # cross-distance matrix for series subset (notice the two-list input)
#' sbD <- proxy::dist(CharTraj[1:10], CharTraj[1:10], method = "SBD", znorm = TRUE)
#'
SBD <- function(x, y, znorm = FALSE, error.check = TRUE) {
    if (is_multivariate(list(x, y))) stop("SBD does not support multivariate series.")

    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }

    nx <- length(x)
    ny <- length(y)

    if (nx > ny) {
        ## The order in which I provide the arguments to NCCc affects 'shift'
        flip <- x
        x <- y
        y <- flip

    } else {
        flip <- NULL
    }

    if (znorm)
        CCseq <- NCCc(zscore(x), zscore(y))
    else
        CCseq <- NCCc(x,y)

    m <- max(CCseq)
    shift <- which.max(CCseq) - max(nx, ny)

    if (is.null(flip)) {
        if (shift < 0L)
            yshift <- y[(-shift + 1L):ny]
        else
            yshift <- c( rep(0, shift), y )

    } else {
        ## Remember, if I flipped them, then I have to shift what is now saved in 'x'
        if (shift < 0L)
            yshift <- c( rep(0, -shift), x )
        else
            yshift <- x[(shift + 1L):ny]
    }

    nys <- length(yshift)

    if (nys < nx)
        yshift <- c( yshift, rep(0, nx-nys) )
    else
        yshift <- yshift[1L:nx]

    list(dist = 1 - m, yshift = yshift)
}

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

SBD_proxy <- function(x, y = NULL, znorm = FALSE, ..., error.check = TRUE, pairwise = FALSE) {
    x <- any2list(x)

    if (error.check) check_consistency(x, "vltslist")
    if (znorm) x <- zscore(x)

    if (is.null(y)) {
        symmetric <- TRUE
        y <- x

        ## Precompute FFTs, padding with zeros as necessary, which will be compensated later
        L <- max(lengths(x)) * 2L - 1L
        fftlen <- stats::nextn(L, 2L)
        fftx <- lapply(x, function(u) { stats::fft(c(u, rep(0, fftlen - length(u)))) })
        ffty <- lapply(fftx, Conj)

    } else {
        symmetric <- FALSE
        y <- any2list(y)
        if (error.check) check_consistency(y, "vltslist")
        if (znorm) y <- zscore(y)

        ## Precompute FFTs, padding with zeros as necessary, which will be compensated later
        L <- max(lengths(x)) + max(lengths(y)) - 1L
        fftlen <- stats::nextn(L, 2L)
        fftx <- lapply(x, function(u) { stats::fft(c(u, rep(0, fftlen - length(u)))) })
        ffty <- lapply(y, function(v) { Conj(stats::fft(c(v, rep(0, fftlen - length(v))))) })
    }

    if (is_multivariate(x) || is_multivariate(y)) stop("SBD does not support multivariate series.")
    pairwise <- isTRUE(pairwise)
    dim_out <- c(length(x), length(y))
    dim_names <- list(names(x), names(y))
    D <- allocate_distmat(length(x), length(y), pairwise, symmetric) ## utils.R

    ## Wrap as needed for foreach
    if (pairwise) {
        x <- split_parallel(x)
        y <- split_parallel(y)
        fftx <- split_parallel(fftx)
        ffty <- split_parallel(ffty)
        validate_pairwise(x, y)
        endpoints <- attr(x, "endpoints")

    } else if (symmetric) {
        endpoints <- symmetric_loop_endpoints(length(x)) ## utils.R
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- x
        fftx <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { fftx })
        ffty <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { ffty })

    } else {
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- split_parallel(y)
        fftx <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { fftx })
        ffty <- split_parallel(ffty)
        endpoints <- attr(y, "endpoints")
    }

    if (bigmemory::is.big.matrix(D)) {
        D_desc <- bigmemory::describe(D)
        noexport <- "D"

    } else {
        D_desc <- NULL
        noexport <- ""
    }

    ## Calculate distance matrix
    foreach(x = x, y = y, fftx = fftx, ffty = ffty, endpoints = endpoints,
            .combine = c,
            .multicombine = TRUE,
            .packages = c("dtwclust", "bigmemory"),
            .export = "sbd_loop",
            .noexport = noexport) %op% {
                bigmat <- !is.null(D_desc)
                d <- if (bigmat) bigmemory::attach.big.matrix(D_desc)@address else D
                sbd_loop(d, x, y, fftx, ffty, fftlen, symmetric, pairwise, endpoints, bigmat)
            }

    D <- D[,]
    if (pairwise) {
        class(D) <- "pairdist"

    } else {
        if (is.null(dim(D))) dim(D) <- dim_out
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }

    attr(D, "method") <- "SBD"
    ## return
    D
}

# ==================================================================================================
# Wrapper for C++
# ==================================================================================================

sbd_loop <- function(d, x, y, fftx, ffty, fftlen, symmetric, pairwise, endpoints, bigmat) {
    .Call(C_sbd_loop,
          d, x, y, fftx, ffty, fftlen, symmetric, pairwise, endpoints, bigmat,
          PACKAGE = "dtwclust")
}
