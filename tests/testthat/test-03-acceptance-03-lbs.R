# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

lb_comp <- function(lower, greater, ...) {
    Map(a = lower, b = greater, MoreArgs = list(...), f = function(a, b, ...) {
        ret <- isTRUE(all.equal(a, b, ...))
        if (!ret) ret <- a < b
        if (!ret) ret <- paste(a, "is NOT less than or equal to", b)
        ret
    })
}

# ==================================================================================================
# L1 norm
# ==================================================================================================

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

    sapply(lb_comp(lbk, lbi), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbk, dtwd), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbi, dtwd), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })

    sapply(lb_comp(lbks, lbis), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbks, dtwds), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbis, dtwds), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })

    expect_equal(lbk, lbks[1L, 1L], info = "single-vs-proxy")
    expect_equal(lbi, lbis[1L, 1L], info = "single-vs-proxy")
    expect_equal(dtwd, dtwds[1L, 1L], info = "single-vs-proxy")
})

# ==================================================================================================
# L2 norm
# ==================================================================================================

test_that("Lower bounds with L2 norm are always leq than DTW.", {
    lbk <- lb_keogh(data_matrix[2L, ], data_matrix[1L, ], norm = "L2", window.size = 15L)$d
    lbi <- lb_improved(data_matrix[2L, ], data_matrix[1L, ], norm = "L2", window.size = 15L)
    dtwd <- dtw_basic(data_matrix[2L, ], data_matrix[1L, ], norm = "L2", window.size = 15L)

    lbks <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                        method = "lbk", window.size = 15L, norm = "L2")
    lbis <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                        method = "lbi", window.size = 15L, norm = "L2")
    dtwds <- proxy::dist(data_reinterpolated[-1L], data_reinterpolated[1L],
                         method = "dtw_basic", window.size = 15L, norm = "L2")

    sapply(lb_comp(lbk, lbi), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbk, dtwd), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbi, dtwd), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })

    sapply(lb_comp(lbks, lbis), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbks, dtwds), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })
    sapply(lb_comp(lbis, dtwds), function(comparison_result) {
        expect_true(comparison_result, info = paste("Result was: ", comparison_result) )
    })

    expect_equal(lbk, lbks[1L, 1L], info = "single-vs-proxy")
    expect_equal(lbi, lbis[1L, 1L], info = "single-vs-proxy")
    expect_equal(dtwd, dtwds[1L, 1L], info = "single-vs-proxy")
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
