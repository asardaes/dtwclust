# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# consistency of dtwb
# ==================================================================================================

test_that("dtw_basic gives the same results as dtw/dtw2", {
    for (rep in 1L:100L) {
        i <- sample(length(CharTraj), 1L)
        j <- sample(length(CharTraj), 1L)
        norm <- sample(c("L1", "L2"), 1L)

        if (i <= 50L)
            window.size <- 20L
        else
            window.size <- NULL

        if (j <= 50L)
            step.pattern <- dtw::symmetric1
        else
            step.pattern <- dtw::symmetric2

        if (is.null(window.size))
            window.type <- "none"
        else
            window.type <- "slantedband"

        if (sample(2L, 1L) == 1L) {
            x <- CharTraj[[i]]
            y <- CharTraj[[j]]
        }
        else {
            x <- CharTrajMV[[i]]
            y <- CharTrajMV[[j]]
        }

        if (norm == "L1")
            suppressWarnings(
                d1 <- dtw(x, y,
                          dist.method = "L1",
                          step.pattern = step.pattern,
                          window.type = window.type,
                          window.size = window.size,
                          keep.internals = TRUE)
            )
        else
            d1 <- dtw2(x, y,
                       step.pattern = step.pattern,
                       window.type = window.type,
                       window.size = window.size,
                       keep.internals = TRUE)

        storage.mode(d1$index1) <- "integer"
        storage.mode(d1$index2) <- "integer"

        d2 <- dtw_basic(x, y, window.size = window.size, norm = norm,
                        step.pattern = step.pattern, backtrack = TRUE)

        expect_equal(d1$distance, d2$distance, info = paste("Indices:", i, j, "- Norm:", norm))

        if (identical(step.pattern, dtw::symmetric2)) {
            d3 <- dtw_basic(x, y, window.size = window.size, norm = norm,
                            step.pattern = step.pattern, normalize = TRUE)

            expect_equal(d1$normalizedDistance, d3, info = paste("Indices:", i, j, "- Norm:", norm))
        }

        expect_identical(d1$index1, d2$index1, info = paste("Indices:", i, j, "- Norm:", norm))
        expect_identical(d1$index2, d2$index2, info = paste("Indices:", i, j, "- Norm:", norm))
    }

    D1_L1 <- proxy::dist(CharTraj[31L:46L], CharTraj[71L:86L], method = "dtw")
    D2_L1 <- proxy::dist(CharTraj[31L:46L], CharTraj[71L:86L], method = "dtw_basic")

    expect_equal(D1_L1, D2_L1, info = "dtw vs dtw_basic", ignore_attr = TRUE)

    D1_L2 <- proxy::dist(CharTraj[31L:46L], CharTraj[71L:86L], method = "dtw2")
    D2_L2 <- proxy::dist(CharTraj[31L:16L], CharTraj[71L:16L], method = "dtw_basic", norm = "L2")

    expect_equal(D1_L1, D2_L1, info = "dtw2 vs dtw_basic", ignore_attr = TRUE)

    D1_L1 <- proxy::dist(CharTrajMV[31L:46L], CharTrajMV[71L:86L], method = "dtw", dist.method = "L1")
    D2_L1 <- proxy::dist(CharTrajMV[31L:46L], CharTrajMV[71L:86L], method = "dtw_basic")

    expect_equal(D1_L1, D2_L1, info = "dtw vs dtw_basic (multivariate)", ignore_attr = TRUE)

    D1_L2 <- proxy::dist(CharTrajMV[31L:46L], CharTrajMV[71L:86L], method = "dtw2")
    D2_L2 <- proxy::dist(CharTrajMV[31L:16L], CharTrajMV[71L:16L], method = "dtw_basic", norm = "L2")

    expect_equal(D1_L1, D2_L1, info = "dtw2 vs dtw_basic (multivariate)", ignore_attr = TRUE)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
