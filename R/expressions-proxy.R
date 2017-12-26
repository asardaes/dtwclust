# ==================================================================================================
# Expressions that are evaluated in the proxy version of some included distances
# ==================================================================================================

# wrap x and y in appropriately arranged lists depending on the kind of distance matrix calculation
foreach_wrap_expression <- expression({
    if (pairwise) {
        x <- split_parallel(x)
        y <- split_parallel(y)
        validate_pairwise(x, y)
        endpoints <- attr(x, "endpoints")

    } else if (symmetric) {
        endpoints <- symmetric_loop_endpoints(length(x)) # utils.R
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- x

    } else {
        x <- lapply(1L:(foreach::getDoParWorkers()), function(dummy) { x })
        y <- split_parallel(y)
        endpoints <- attr(y, "endpoints")
    }
})

# perform the actual calculation of the distance matrix with support for parallelization
#' @importFrom bigmemory attach.big.matrix
#' @importFrom bigmemory is.big.matrix
#'
foreach_loop_expression <- expression({
    export <- c("enlist")

    do.call(what = foreach::foreach, quote = TRUE, args = enlist(
        x = x, y = y, endpoints = endpoints, dots = foreach_extra_args, # specified by some
        .combine = c, .multicombine = TRUE, .packages = packages,
        .export = export, .noexport = noexport
    )) %op% {
        dots <- lapply(dots, eval, envir = environment())
        bigmat <- !is.null(D_desc)
        d <- if (bigmat) bigmemory::attach.big.matrix(D_desc)@address else D
        do.call(.distfun_, # specified by each function
                enlist(d = d,
                       x = x,
                       y = y,
                       symmetric = symmetric,
                       pairwise = pairwise,
                       endpoints = endpoints,
                       bigmat = bigmat,
                       dots = dots),
                TRUE)
    }

    if (pairwise || bigmemory::is.big.matrix(D))
        D <- D[,] # coerce: column vector -> vector || big.matrix -> R matrix
})
