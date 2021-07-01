context("    Miscellaneous functions")

# ==================================================================================================
# setup
# ==================================================================================================

# Original objects in env
ols <- ls()

test_that("Corner cases for check_consistency.", {
    expect_error(dtwclust:::check_consistency(list(), "tslist"), "empty")
    expect_error(dtwclust:::check_consistency(list(), "vltslist"), "empty")
    expect_error(dtwclust:::check_consistency(NULL, "window"), "provide")
    expect_error(dtwclust:::check_consistency(-1L, "window"), "negative")
    expect_warning(
        dtwclust:::check_consistency("foo", "cent",
                                     clus_type = "hierarchical",
                                     cent_char = "bar",
                                     cent_missing = FALSE),
        "ignored"
    )
})

# ==================================================================================================
# compute_envelope
# ==================================================================================================

test_that("compute_envelope catches input errors.", {
    expect_error(compute_envelope(matrix(0L,2L,2L)), "univariate")
    expect_error(compute_envelope(1:3, 7L), "Window")
})

# ==================================================================================================
# reinterpolate
# ==================================================================================================

test_that("reinterpolate function works correctly for supported inputs.", {
    expect_length(reinterpolate(data[[1L]], 200L), 200L)
    expect_length(reinterpolate(data[[1L]], 200L, TRUE), 200L)

    expect_identical(ncol(reinterpolate(data_matrix, 200L, multivariate = FALSE)), 200L)
    expect_identical(nrow(reinterpolate(data_matrix, 200L, multivariate = TRUE)), 200L)

    data_frame <- base::as.data.frame(data_matrix)
    expect_identical(ncol(reinterpolate(data_frame, 200L, multivariate = FALSE)), 200L)
    expect_identical(nrow(reinterpolate(data_frame, 200L, multivariate = TRUE)), 200L)
})

# ==================================================================================================
# NCCc
# ==================================================================================================

test_that("NCCc function works correctly for supported inputs.", {
    expect_true(is.infinite(NCCc(0, 1:2)))

    x_uv <- data[[1L]]
    x_mv <- data_multivariate[[1L]]
    invalid_inputs <- list(
        "NA" = NA,
        "NULL" = NULL,
        "char" = "1",
        "bool" = TRUE,
        "empty" = numeric(),
        "miss" = c(1, NA)
    )

    expect_error(NCCc(x_uv, x_mv), regexp = "dimension")
    expect_error(NCCc(x_mv, x_uv), regexp = "dimension")
    expect_error(NCCc(x_mv, x_mv), regexp = "multivariate")
    for (input in names(invalid_inputs)) {
        expect_error(do.call("NCCc", list(x = invalid_inputs[[input]], y = x_uv)),
                     info = paste("with", input, "input in x"))
        expect_error(do.call("NCCc", list(x = x_uv, y = invalid_inputs[[input]])),
                     info = paste("with", input, "input in y"))
    }
    expect_length(NCCc(x_uv, x_uv), 2L * length(x_uv) - 1L)
})

# ==================================================================================================
# zscore
# ==================================================================================================

test_that("zscore function works correctly for supported inputs.", {
    normalized_series <- zscore(data[[1L]])
    expect_length(normalized_series, length(data[[1L]]))
    expect_true(!is.null(attributes(zscore(data[[1L]], keep.attributes = TRUE))))

    normalized_matrix <- zscore(data_matrix, multivariate = FALSE, keep.attributes = TRUE)
    expect_identical(ncol(normalized_matrix), ncol(data_matrix))
    expect_true(!is.null(attributes(normalized_matrix)))

    normalized_matrix <- zscore(data_matrix, multivariate = TRUE, keep.attributes = TRUE)
    expect_identical(nrow(normalized_matrix), nrow(data_matrix))
    expect_true(!is.null(attributes(normalized_matrix)))

    data_frame <- base::as.data.frame(data_matrix)

    normalized_data_frame <- zscore(data_frame, multivariate = FALSE, keep.attributes = TRUE)
    expect_identical(ncol(normalized_data_frame), ncol(data_frame))
    expect_true(!is.null(attributes(normalized_data_frame)))

    normalized_data_frame <- zscore(data_frame, multivariate = TRUE, keep.attributes = TRUE)
    expect_identical(nrow(normalized_data_frame), nrow(data_frame))
    expect_true(!is.null(attributes(normalized_data_frame)))
})

# ==================================================================================================
# %op%
# ==================================================================================================

test_that("%op% catches errors as expected.", {
    skip_if(foreach::getDoParWorkers() %% 2L == 1L, "sequential case")
    prev_threads <- foreach(dummy = 1L:foreach::getDoParWorkers()) %dopar% {
        nthreads <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "1"))
        RcppParallel::setThreadOptions()
        nthreads
    }
    expect_error(
        dtwclust:::`%op%`(foreach(i = 1L:2L, .packages = "dtwclust"), {
            if (i == 1L) stop("test") else NULL
        })
    )
    workers_threads <- foreach(dummy = 1L:foreach::getDoParWorkers()) %dopar% {
        Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
    }
    sapply(workers_threads, function(nthreads_env_var) {
        expect_false(nzchar(nthreads_env_var))
    })
    foreach(nthreads = prev_threads) %dopar% {
        RcppParallel::setThreadOptions(nthreads)
        NULL
    }
})

# ==================================================================================================
# PairTracker
# ==================================================================================================

test_that("PairTracker works as expected.", {
    expect_error(dtwclust:::PairTracker$new(0L), "size")

    set.seed(103289L)
    tracker <- dtwclust:::PairTracker$new(3L)
    expect_error(tracker$link(0L, 1L, 0L), "Invalid indices")
    expect_error(tracker$link(1L, 2L, 2L), "link type")

    expect_false(tracker$link(1L, 2L, -1L), "link 1 and 2 as dont-know")
    invisible(sapply(1L:100L, function(dummy) {
        expect_false(all(c(1L,2L) %in% tracker$get_unseen_pair()),
                     "1 and 2 are linked, so they should not appear again")
        NULL
    }))

    expect_false(tracker$link(2L, 3L, 0L), "link 2 and 3 as cannot-link")
    invisible(sapply(1L:100L, function(dummy) {
        expect_true(all(c(1L,3L) %in% tracker$get_unseen_pair()),
                    "1 and 3 are the only unlinked indices")
        NULL
    }))

    expect_false(tracker$link(1L, 3L, 1L), "link 1 and 3 as must-link")
    expect_null(tracker$get_unseen_pair(), "aggregate graph should be complete now")
    expect_true(tracker$link(2L, 3L, 1L), "must-link graph should be connected now")
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
