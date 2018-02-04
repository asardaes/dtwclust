#' @importFrom parallel nextRNGStream
#'
rng_seq <- function(n, seed = NULL, simplify = FALSE) {
    stopifnot(RNGkind()[1L] == rng_kind, n > 0L)

    if (!is.null(seed)) {
        prev_seed <- get(".Random.seed", globalenv())
        if (length(seed) == 1L) {
            set.seed(as.integer(seed))
        }
        else if (length(seed) == 7L) {
            assign(".Random.seed", seed, globalenv())
        }
        else {
            stop("Invalid seed provided") # nocov
        }
        on.exit(assign(".Random.seed", prev_seed, globalenv()))
    }

    first_seed <- get(".Random.seed", globalenv())
    if (n == 1L) {
        if (!isTRUE(simplify)) first_seed <- list(first_seed)
        return(first_seed) # return 1
    }

    seed_seq <- vector("list", n)
    seed_seq[[1L]] <- first_seed
    for (i in 2L:n) seed_seq[[i]] <- parallel::nextRNGStream(seed_seq[[i-1L]])
    # return
    seed_seq
}
