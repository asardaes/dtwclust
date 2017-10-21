context("\tMiscellaneous functions")

# ==================================================================================================
# setup
# ==================================================================================================

# Original objects in env
ols <- ls()

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
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
