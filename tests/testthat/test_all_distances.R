context("Test included distances")

ctrl <- new("dtwclustControl", window.size = 18L)
x <- data_reinterpolated[1L:20L]
centroids <- x[c(1L, 15L)]
attr(centroids, "id_cent") <- c(1L, 15L)

# =================================================================================================
# distance functions
# =================================================================================================

test_that("Operations with dtwclustFamily@dist give expected results", {
     ## ---------------------------------------------------------- L2
     distmat <- proxy::dist(x, x, method = "L2")

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "L2")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

     expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                  tolerance = 0, check.attributes = FALSE)

     class(sub_distmat) <- "matrix"

     expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                  tolerance = 0, check.attributes = FALSE)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   distmat = distmat,
                   dist = "L2")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

     expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                  tolerance = 0, check.attributes = FALSE)

     expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                  tolerance = 0, check.attributes = FALSE)

     ## ---------------------------------------------------------- lbk
     distmat <- proxy::dist(x, x, method = "lbk", window.size = 18L)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "lbk")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

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

     ## ---------------------------------------------------------- lbi
     distmat <- proxy::dist(x, x, method = "lbi", window.size = 18L)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "lbi")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

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

     ## ---------------------------------------------------------- sbd
     distmat <- proxy::dist(x, x, method = "sbd")

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "sbd")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

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

     ## ---------------------------------------------------------- dtw_lb
     distmat <- proxy::dist(x, x, method = "dtw_lb", window.size = 18L)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "dtw_lb")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

     expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                  tolerance = 0, check.attributes = FALSE)

     expect_equal(sub_distmat, proxy::dist(x, centroids, method = "dtw_lb", window.size = 18L),
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

     ## ---------------------------------------------------------- dtw
     distmat <- proxy::dist(x, x, method = "dtw", window.type = "slantedband", window.size = 18L)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "dtw")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

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

     ## ---------------------------------------------------------- dtw2
     distmat <- proxy::dist(x, x, method = "dtw2", window.type = "slantedband", window.size = 18L)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "dtw2")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

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

     ## ---------------------------------------------------------- dtw_basic
     distmat <- proxy::dist(x, x, method = "dtw_basic", window.size = 18L)

     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "dtw_basic")

     whole_distmat <- family@dist(x)
     sub_distmat <- family@dist(x, centroids)

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
})

# =================================================================================================
# invalid distance
# =================================================================================================

test_that("Errors in centroid argument are correctly detected.", {
     expect_error(dtwclust(data_matrix, k = 20, distance = mean), "proxy", info = "Function")

     expect_error(dtwclust(data_matrix, k = 20, distance = "dummy"), "proxy", info = "Unregistered")

     expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                  "different length", info = "LBK")

     expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                  "different length", info = "LBI")

     expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                  "different length", info = "DTW_LB")
})
