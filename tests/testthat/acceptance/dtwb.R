context("\tConsistency of dtw_basic")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# consistency of dtwb
# =================================================================================================

test_that("dtw_basic gives the same results as dtw/dtw2", {
    for (rep in 1L:100L) {
        i <- sample(length(data), 1L)
        j <- sample(length(data), 1L)
        norm <- sample(c("L1", "L2"), 1L)

        if (i <= 50L)
            window.size <- 20L
        else
            window.size <- NULL

        if (j <= 50L)
            step.pattern <- symmetric1
        else
            step.pattern <- symmetric2

        if (is.null(window.size))
            window.type <- "none"
        else
            window.type <- "slantedband"

        x <- data[[i]]
        y <- data[[j]]

        if (norm == "L1")
            d1 <- dtw(x, y, step.pattern = step.pattern,
                      window.type = window.type, window.size = window.size,
                      keep.internals = TRUE)
        else
            d1 <- dtw2(x, y, step.pattern = step.pattern,
                       window.type = window.type, window.size = window.size,
                       keep.internals = TRUE)

        storage.mode(d1$index1) <- "integer"
        storage.mode(d1$index2) <- "integer"

        d2 <- dtw_basic(x, y, window.size = window.size, norm = norm,
                        step.pattern = step.pattern, backtrack = TRUE)

        expect_identical(d1$distance, d2$distance, info = paste("Indices:", i, j, "- Norm:", norm))

        if (identical(step.pattern, symmetric2)) {
            d3 <- dtw_basic(x, y, window.size = window.size, norm = norm,
                            step.pattern = step.pattern, normalize = TRUE)

            expect_identical(d1$normalizedDistance, d3, info = paste("Indices:", i, j, "- Norm:", norm))
        }

        expect_identical(d1$index1, d2$index1, info = paste("Indices:", i, j, "- Norm:", norm))
        expect_identical(d1$index2, d2$index2, info = paste("Indices:", i, j, "- Norm:", norm))
    }

    D1_L1 <- proxy::dist(data[31L:46L], data[71L:86L], method = "dtw")
    D2_L1 <- proxy::dist(data[31L:46L], data[71L:86L], method = "dtw_basic")

    expect_equal(D1_L1, D2_L1, tolerance = 0, check.attributes = FALSE,
                 info = "dtw vs dtw_basic")

    D1_L2 <- proxy::dist(data[31L:46L], data[71L:86L], method = "dtw2")
    D2_L2 <- proxy::dist(data[31L:16L], data[71L:16L], method = "dtw_basic", norm = "L2")

    expect_equal(D1_L1, D2_L1, tolerance = 0, check.attributes = FALSE,
                 info = "dtw2 vs dtw_basic")
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
