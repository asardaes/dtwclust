context("Test included distances")

ctrl <- new("dtwclustControl", window.size = 18L)
x <- data_reinterpolated[1L:20L]
centroids <- x[c(1L, 15L)]
attr(centroids, "id_cent") <- c(1L, 15L)

# =================================================================================================
# distance functions
# =================================================================================================

test_that("Operations with dtwclustFamily@dist give expected results", {
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

     distmat_lbk <- whole_distmat

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

     distmat_lbi <- whole_distmat

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

     distmat_sbd <- whole_distmat

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

     distmat_dtwlb <- whole_distmat

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

     distmat_dtw <- whole_distmat

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

     distmat_dtw2 <- whole_distmat

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

     distmat_dtwb <- whole_distmat

     skip_on_cran()

     expect_equal_to_reference(distmat_lbk, file_name(distmat_lbk), info = "LBK")
     expect_equal_to_reference(distmat_lbi, file_name(distmat_lbi), info = "LBI")
     expect_equal_to_reference(distmat_sbd, file_name(distmat_sbd), info = "SBD")
     expect_equal_to_reference(distmat_dtwlb, file_name(distmat_dtwlb), info = "DTW_LB")
     expect_equal_to_reference(distmat_dtw, file_name(distmat_dtw), info = "DTW")
     expect_equal_to_reference(distmat_dtw2, file_name(distmat_dtw2), info = "DTW2")
     expect_equal_to_reference(distmat_dtwb, file_name(distmat_dtwb), info = "DTW_BASIC")
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
