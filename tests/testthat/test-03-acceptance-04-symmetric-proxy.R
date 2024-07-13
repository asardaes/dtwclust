# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# split_parallel_symmetric in UTILS-utils.R
# ==================================================================================================

test_that("The indices for load balancing always cover the whole distance matrix.", {
    ## random test
    set.seed(as.integer(Sys.time()))

    for (i in seq_len(30L)) {
        num_series <- round(runif(1L, 5, 100))
        num_workers <- round(runif(1L, 3, 100))
        info <- paste("#series = ", num_series, ", #workers = ", num_workers, collapse = "")
        x <- data_reinterpolated[1L:num_series]
        distmat <- matrix(0, num_series, num_series)
        indices <- dtwclust:::split_parallel_symmetric(num_series, num_workers)

        for (ids in indices) {
            if (isTRUE(attr(ids, "trimat"))) {
                ## assign upper part of lower triangular
                ul <- ids$ul
                if (length(ul) > 1L)
                    distmat[ul,ul] <- base::as.matrix(do.call(
                        proxy::dist,
                        list(x = x[ul],
                             y = NULL,
                             method = "L2")
                    ))
                ## assign lower part of lower triangular
                ll <- ids$ll
                if (length(ll) > 1L)
                    distmat[ll,ll] <- base::as.matrix(do.call(
                        proxy::dist,
                        list(x = x[ll],
                             y = NULL,
                             method = "L2")
                    ))
            } else {
                rows <- attr(ids, "rows")
                mat_chunk <- base::as.matrix(do.call(
                    proxy::dist,
                    list(x = x[rows],
                         y = x[ids],
                         method = "L2")
                ))
                ## assign matrix chunks
                distmat[rows,ids] <- mat_chunk
                distmat[ids,rows] <- t(mat_chunk)
            }
        }

        diagonal_is_zero <- all.equal(diag(distmat), rep(0, num_series))
        lt <- all(distmat[lower.tri(distmat)] > 0)
        ut <- all(distmat[upper.tri(distmat)] > 0)
        nondiagonal_is_nonzero <- lt & ut

        expect_true(diagonal_is_zero, info = info)
        expect_true(nondiagonal_is_nonzero, info = info)
    }
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
