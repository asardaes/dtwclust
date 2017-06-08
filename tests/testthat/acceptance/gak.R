context("\tSymmetric GAK")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# symmetric gak
# =================================================================================================

test_that("Symmetric univariate GAK distance gives expected results.", {
    D1 <- proxy::dist(data_subset, method = "gak", window.size = 18L)
    sigma <- attr(D1, "sigma")
    D2 <- proxy::dist(data_subset, data_subset, method = "gak", window.size = 18L, sigma = sigma)
    D3 <- proxy::dist(data_subset, data_subset, method = "gak",window.size = 18L, sigma = sigma, symmetric = TRUE)
    D4 <- sapply(data_subset, GAK, y = data_subset[[1L]], window.size = 18L, sigma = sigma)

    expect_equal(D1, D2, info = "single-double-arg", tolerance = 0, check.attributes = FALSE)
    expect_equal(D2, D3, info = "forced-symmetry", tolerance = 0, check.attributes = FALSE)
    expect_equal(D1[ , 1L], D4, info = "manual-vs-proxy", tolerance = 0, check.attributes = FALSE)
})

test_that("Symmetric multivariate GAK distance gives expected results.", {
    D1 <- proxy::dist(data_multivariate, method = "gak", window.size = 18L)
    sigma <- attr(D1, "sigma")
    D2 <- proxy::dist(data_multivariate, data_multivariate, method = "gak", window.size = 18L, sigma = sigma)
    D3 <- proxy::dist(data_multivariate, data_multivariate, method = "gak", window.size = 18L, sigma = sigma, symmetric = TRUE)
    D4 <- sapply(data_multivariate, GAK, y = data_multivariate[[1L]], window.size = 18L, sigma = sigma)

    expect_equal(D1, D2, info = "single-double-arg", tolerance = 0, check.attributes = FALSE)
    expect_equal(D2, D3, info = "forced-symmetry", tolerance = 0, check.attributes = FALSE)
    expect_equal(D1[ , 1L], D4, info = "manual-vs-proxy", tolerance = 0, check.attributes = FALSE)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
