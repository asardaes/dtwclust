context("    Proxy distances")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

## Data
x <- data_reinterpolated[3L:8L]

# ==================================================================================================
# proxy distances
# ==================================================================================================

test_that("Included proxy distances can be called and give expected dimensions.", {
    for (distance in dtwclust:::distances_included) {
        d <- proxy::dist(x, method = distance, window.size = 15L, sigma = 100, normalize = TRUE)
        expect_identical(dim(d), c(length(x), length(x)), info = paste(distance, "single-arg"))

        d2 <- proxy::dist(x, x, method = distance, window.size = 15L, sigma = 100, normalize = TRUE)
        expect_equal(d2, d, check.attributes = FALSE,
                     info = paste(distance, "double-arg"))

        d3 <- proxy::dist(x[1L], x, method = distance, window.size = 15L, sigma = 100, normalize = TRUE)
        class(d3) <- c("matrix", "array")
        expect_identical(dim(d3), c(1L, length(x)), info = paste(distance, "one-vs-many"))

        d4 <- proxy::dist(x, x[1L], method = distance, window.size = 15L, sigma = 100, normalize = TRUE)
        class(d4) <- c("matrix", "array")
        expect_identical(dim(d4), c(length(x), 1L), info = paste(distance, "many-vs-one"))

        # dtw_lb will give different results below because of how it works
        if (distance == "dtw_lb") next

        expect_equal(d3, d[1L, , drop = FALSE], check.attributes = FALSE,
                     info = paste(distance, "one-vs-many-vs-distmat"))
        expect_equal(d4, d[ , 1L, drop = FALSE], check.attributes = FALSE,
                     info = paste(distance, "many-vs-one-vs-distmat"))

        dots <- list()
        if (distance %in% c("lb_keogh", "lb_improved"))
            dots <- list(window.size = 15L)
        else if (distance %in% c("gak"))
            dots <- list(window.size = 15L, sigma = 100)
        else if (distance %in% c("dtw_basic"))
            dots <- list(window.size = 15L, normalize = TRUE)

        manual_distmat <- sapply(x, function(j) {
            sapply(x, function(i) {
                d <- do.call(distance, dtwclust:::enlist(x = i, y = j, dots = dots), TRUE)
                if (distance %in% c("lb_keogh", "sbd")) d <- d$d
                d
            })
        })
        expect_equal(as.matrix(d), manual_distmat, check.attributes = FALSE,
                     info = paste("manual distmat vs proxy version using", distance))
    }
})

test_that("Parameter errors in included distances are detected.", {
    expect_error(proxy::dist(data_multivariate, method = "dtw_lb"), "multivariate")
    expect_error(proxy::dist(list(), method = "dtw_lb"), "Empty")
    expect_error(proxy::dist(data_subset, list(), method = "dtw_lb"), "Empty")
    expect_error(proxy::dist(data_subset, method = "sdtw", gamma = -1))
    expect_error(proxy::dist(data_subset, method = "gak", sigma = -1))
    expect_error(proxy::dist(data_subset, method = "dtw_basic", step.pattern = dtw::asymmetric))
    expect_error(proxy::dist(data_subset, method = "dtw_basic",
                             step.pattern = dtw::symmetric1, normalize = TRUE))
})

# ==================================================================================================
# proxy pairwise distances
# ==================================================================================================

test_that("Included proxy distances can be called for pairwise = TRUE and give expected length", {
    for (distance in dtwclust:::distances_included) {
        ## sbd doesn't always return zero, so tolerance is left alone here

        d <- proxy::dist(x, method = distance,
                         window.size = 15L, step.pattern = dtw::symmetric1,
                         pairwise = TRUE)
        class(d) <- "numeric"
        expect_null(dim(d), paste("distance =", distance))
        expect_identical(length(d), length(x), info = paste(distance, "pairwise single-arg"))
        if (distance != "sdtw")
            expect_equal(d, rep(0, length(d)), check.attributes = FALSE,
                         info = paste(distance, "pairwise single all zero"))

        d2 <- proxy::dist(x, x, method = distance,
                          window.size = 15L, step.pattern = dtw::symmetric1,
                          pairwise = TRUE)
        class(d2) <- "numeric"
        expect_null(dim(d2), paste("distance =", distance))
        expect_identical(length(d2), length(x), info = paste(distance, "pairwise double-arg"))
        if (distance != "sdtw")
            expect_equal(d, rep(0, length(d2)), check.attributes = FALSE,
                         info = paste(distance, "pairwise double all zero"))

        expect_error(proxy::dist(x[1L:3L], x[4L:5L], method = distance,
                                 window.size = 15L, pairwise = TRUE),
                     "same amount",
                     info = paste(distance, "invalid pairwise"))
    }
})

# ==================================================================================================
# proxy similarities
# ==================================================================================================

test_that("Included proxy similarities can be called and give expected dimensions.", {
    for (distance in c("uGAK")) {
        d <- proxy::simil(x, method = distance, sigma = 100)
        expect_identical(dim(d), c(length(x), length(x)), info = paste(distance, "single-arg"))

        d2 <- proxy::simil(x, x, method = distance, sigma = 100)
        expect_equal(d2, d, check.attributes = FALSE,
                     info = paste(distance, "double-arg"))

        d3 <- proxy::simil(x[1L], x, method = distance, sigma = 100)
        class(d3) <- c("matrix", "array")
        expect_identical(dim(d3), c(1L, length(x)), info = paste(distance, "one-vs-many"))

        d4 <- proxy::simil(x, x[1L], method = distance, sigma = 100)
        class(d4) <- c("matrix", "array")
        expect_identical(dim(d4), c(length(x), 1L), info = paste(distance, "many-vs-one"))

        expect_equal(d3, d[1L, , drop = FALSE], check.attributes = FALSE,
                     info = paste(distance, "one-vs-many-vs-distmat"))
        expect_equal(d4, d[ , 1L, drop = FALSE], check.attributes = FALSE,
                     info = paste(distance, "many-vs-one-vs-distmat"))
    }
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
