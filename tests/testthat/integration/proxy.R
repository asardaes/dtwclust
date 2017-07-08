context("\tProxy distances")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

## Data
x <- data_reinterpolated[3L:8L]

# =================================================================================================
# proxy distances
# =================================================================================================

test_that("Included proxy distances can be called and give expected dimensions.", {
    for (distance in dtwclust:::distances_included) {
        d <- proxy::dist(x, method = distance, window.size = 15L, sigma = 100)
        expect_identical(dim(d), c(length(x), length(x)), info = paste(distance, "single-arg"))

        d2 <- proxy::dist(x, x, method = distance, window.size = 15L, sigma = 100)
        expect_equal(d2, d, check.attributes = FALSE,
                     info = paste(distance, "double-arg"))

        d3 <- proxy::dist(x[1L], x, method = distance, window.size = 15L, sigma = 100)
        class(d3) <- "matrix"
        expect_identical(dim(d3), c(1L, length(x)), info = paste(distance, "one-vs-many"))

        d4 <- proxy::dist(x, x[1L], method = distance, window.size = 15L, sigma = 100)
        class(d4) <- "matrix"
        expect_identical(dim(d4), c(length(x), 1L), info = paste(distance, "many-vs-one"))

        ## dtw_lb will give different results below because of how it works
        if (distance == "dtw_lb") next

        expect_equal(d3, d[1L, , drop = FALSE], check.attributes = FALSE,
                     info = paste(distance, "one-vs-many-vs-distmat"))
        expect_equal(d4, d[ , 1L, drop = FALSE], check.attributes = FALSE,
                     info = paste(distance, "many-vs-one-vs-distmat"))
    }
})

# =================================================================================================
# proxy pairwise distances
# =================================================================================================

test_that("Included proxy distances can be called for pairwise = TRUE and give expected length", {
    for (distance in dtwclust:::distances_included) {
        ## sbd doesn't always return zero, so tolerance is left alone here

        d <- proxy::dist(x, method = distance, window.size = 15L, pairwise = TRUE)
        class(d) <- "numeric"
        expect_null(dim(d))
        expect_identical(length(d), length(x), info = paste(distance, "pairwise single-arg"))
        expect_equal(d, rep(0, length(d)), check.attributes = FALSE,
                     info = paste(distance, "pairwise single all zero"))

        d2 <- proxy::dist(x, x, method = distance, window.size = 15L, pairwise = TRUE)
        class(d2) <- "numeric"
        expect_null(dim(d2))
        expect_identical(length(d2), length(x), info = paste(distance, "pairwise double-arg"))
        expect_equal(d, rep(0, length(d2)), check.attributes = FALSE,
                     info = paste(distance, "pairwise double all zero"))

        expect_error(proxy::dist(x[1L:3L], x[4L:5L], method = distance,
                                 window.size = 15L, pairwise = TRUE),
                     "same amount",
                     info = paste(distance, "invalid pairwise"))
    }
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
