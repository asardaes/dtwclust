context("\tConsistency of DTW lower bounds")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# L1 norm
# =================================================================================================

test_that("Lower bounds with L1 norm are always leq than DTW.", {
    lbk <- lb_keogh(data_matrix[2L, ], data_matrix[1L, ], window.size = 15L)$d
    lbi <- lb_improved(data_matrix[2L, ], data_matrix[1L, ], window.size = 15L)
    dtwd <- dtw_basic(data_matrix[2L, ], data_matrix[1L, ], window.size = 15L)

    lbks <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                        method = "lbk", window.size = 15L)
    lbis <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                        method = "lbi", window.size = 15L)
    dtwds <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                         method = "dtw_basic", window.size = 15L)

    expect_true(lbk <= lbi, info = "single")
    expect_true(lbk <= dtwd, info = "single")
    expect_true(lbi <= dtwd, info = "single")

    expect_true(all(lbks <= lbis), info = "multiple")
    expect_true(all(lbks <= dtwds), info = "multiple")
    expect_true(all(lbis <= dtwds), info = "multiple")

    expect_true(lbk == lbks[1L, 1L], info = "single-vs-proxy")
    expect_true(lbi == lbis[1L, 1L], info = "single-vs-proxy")
    expect_true(dtwd == dtwds[1L, 1L], info = "single-vs-proxy")
})

# =================================================================================================
# L2 norm
# =================================================================================================

test_that("Lower bounds with L2 norm are always leq than DTW.", {
    lbk <- lb_keogh(data_matrix[2L, ], data_matrix[1L, ], norm = "L2", window.size = 15L)$d
    lbi <- lb_improved(data_matrix[2L, ], data_matrix[1L, ], norm = "L2", window.size = 15L)
    dtwd <- dtw_basic(data_matrix[2L, ], data_matrix[1L, ], norm = "L2", window.size = 15L)

    lbks <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                        method = "lbk", window.size = 15L,
                        norm = "L2")
    lbis <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                        method = "lbi", window.size = 15L,
                        norm = "L2")
    dtwds <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                         method = "dtw_basic", window.size = 15L,
                         norm = "L2")

    expect_true(lbk <= lbi, info = "single")
    expect_true(lbk <= dtwd, info = "single")
    expect_true(lbi <= dtwd, info = "single")

    expect_true(all(lbks <= lbis), info = "multiple")
    expect_true(all(lbks <= dtwds), info = "multiple")
    expect_true(all(lbis <= dtwds), info = "multiple")

    expect_true(lbk == lbks[1L, 1L], info = "single-vs-proxy")
    expect_true(lbi == lbis[1L, 1L], info = "single-vs-proxy")
    expect_true(dtwd == dtwds[1L, 1L], info = "single-vs-proxy")
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
