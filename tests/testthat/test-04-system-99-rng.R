test_that("The RNGkind was not affected by dtwclust.", {
    expect_identical(RNGkind()[1L], default_rngkind)
})
