# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# symmetric gak
# ==================================================================================================

test_that("Symmetric univariate GAK distance gives expected results.", {
    D1 <- proxy::dist(data_subset, method = "gak", window.size = 18L)
    sigma <- attr(D1, "sigma")
    D2 <- proxy::dist(data_subset, data_subset, method = "gak", window.size = 18L, sigma = sigma)
    D3 <- sapply(data_subset, GAK, y = data_subset[[1L]], window.size = 18L, sigma = sigma)

    expect_equal(as.matrix(D1), unclass(D2), info = "single-double-arg", ignore_attr = TRUE)
    expect_equal(c(0, D1[1L:(length(D3) - 1L)]), D3, info = "manual-vs-proxy", ignore_attr = TRUE)
})

test_that("Symmetric multivariate GAK distance gives expected results.", {
    D1 <- proxy::dist(data_multivariate, method = "gak", window.size = 18L)
    sigma <- attr(D1, "sigma")
    D2 <- proxy::dist(data_multivariate, data_multivariate, method = "gak", window.size = 18L, sigma = sigma)
    D3 <- sapply(data_multivariate, GAK, y = data_multivariate[[1L]], window.size = 18L, sigma = sigma)

    expect_equal(as.matrix(D1), unclass(D2), info = "single-double-arg", ignore_attr = TRUE)
    expect_equal(c(0, D1[1L:(length(D3) - 1L)]), D3, info = "manual-vs-proxy", ignore_attr = TRUE)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
