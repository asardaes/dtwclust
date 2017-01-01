context("\tFamilies and proxy distances")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

ctrl <- new("dtwclustControl", window.size = 18L)
x <- data_reinterpolated_subset
centroids <- x[c(1L, 15L)]
attr(centroids, "id_cent") <- c(1L, 15L)

# =================================================================================================
# distance functions
# =================================================================================================

test_that("Operations with dtwclustFamily@dist give expected results", {
    ## ---------------------------------------------------------- lbk
    distmat <- proxy::dist(x, x, method = "lbk", window.size = ctrl@window.size)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "lbk")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "lbk")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_lbk <- whole_distmat

    ## ---------------------------------------------------------- lbi
    distmat <- proxy::dist(x, x, method = "lbi", window.size = ctrl@window.size)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "lbi")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "lbi")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_lbi <- whole_distmat

    ## ---------------------------------------------------------- sbd
    distmat <- proxy::dist(x, x, method = "sbd")

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "sbd")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "sbd")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_sbd <- whole_distmat

    ## ---------------------------------------------------------- dtw_lb
    distmat <- proxy::dist(x, x, method = "dtw_lb", window.size = ctrl@window.size)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "dtw_lb")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, proxy::dist(x, centroids, method = "dtw_lb", window.size = ctrl@window.size),
                 info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "dtw_lb")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_dtwlb <- whole_distmat

    ## ---------------------------------------------------------- dtw
    distmat <- proxy::dist(x, x, method = "dtw", window.type = "slantedband", window.size = ctrl@window.size)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "dtw")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "dtw")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_dtw <- whole_distmat

    ## ---------------------------------------------------------- dtw2
    distmat <- proxy::dist(x, x, method = "dtw2", window.type = "slantedband", window.size = ctrl@window.size)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "dtw2")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "dtw2")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_dtw2 <- whole_distmat

    ## ---------------------------------------------------------- dtw_basic
    distmat <- proxy::dist(x, x, method = "dtw_basic", window.size = ctrl@window.size)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "dtw_basic")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "dtw_basic")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_dtwb <- whole_distmat

    ## ---------------------------------------------------------- gak
    distmat <- proxy::dist(x, x, method = "gak",
                           window.size = ctrl@window.size,
                           sigma = 100)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "gak")

    whole_distmat <- family@dist(x, sigma = 100)
    sub_distmat <- family@dist(x, centroids, sigma = 100)
    pdist <- family@dist(x, sigma = 100, pairwise = TRUE)

    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    class(sub_distmat) <- "matrix"

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    family <- new("dtwclustFamily",
                  control = ctrl,
                  distmat = distmat,
                  dist = "gak")

    whole_distmat <- family@dist(x, sigma = 100)
    sub_distmat <- family@dist(x, centroids, sigma = 100)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    distmat_gak <- whole_distmat

    ## ---------------------------------------------------------- for references
    assign("distmat_lbk", distmat_lbk, persistent)
    assign("distmat_lbi", distmat_lbi, persistent)
    assign("distmat_sbd", distmat_sbd, persistent)
    assign("distmat_dtwlb", distmat_dtwlb, persistent)
    assign("distmat_dtw", distmat_dtw, persistent)
    assign("distmat_dtw2", distmat_dtw2, persistent)
    assign("distmat_dtwb", distmat_dtwb, persistent)
    assign("distmat_gak", distmat_gak, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
