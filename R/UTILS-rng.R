#' @importFrom parallel nextRNGStream
#'
rng_seq <- function(n, seed = NULL, simplify = FALSE) {
    stopifnot(RNGkind()[1L] == dtwclust_rngkind, n > 0L)

    if (!is.null(seed)) {
        prev_seed <- get(".Random.seed", .GlobalEnv)
        if (length(seed) == 1L) {
            set.seed(as.integer(seed))
        }
        else if (length(seed) == 7L) {
            assign(".Random.seed", seed, .GlobalEnv)
        }
        else {
            stop("Invalid seed provided") # nocov
        }
        on.exit(assign(".Random.seed", prev_seed, .GlobalEnv))
    }

    first_seed <- get(".Random.seed", .GlobalEnv)
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

# see https://stackoverflow.com/a/20998531/5793905 and the link there
handle_rngkind <- function() {
    current_seed <- get(".Random.seed", .GlobalEnv)
    previous_rngkind <- RNGkind(dtwclust_rngkind)[1L]
    if (previous_rngkind != dtwclust_rngkind) {
        # evaluate on.exit on the caller's environment, otherwise it would execute immediately
        do.call(on.exit, envir = parent.frame(), args = list(substitute({
            RNGkind(previous_rngkind)
            assign(".Random.seed", current_seed, .GlobalEnv)
        })))
    }
    invisible()
}
