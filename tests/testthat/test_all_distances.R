context("Test included distances")

# =================================================================================================
# proxy distances
# =================================================================================================

included_distances <- c("lbk", "lbi", "sbd", "dtw_basic", "dtw_lb")
x <- data_reinterpolated[3L:8L]

test_that("Included proxy distances can be called and give expected dimensions.", {
     for (distance in included_distances) {
          d <- proxy::dist(x, method = distance, window.size = 15L)

          expect_identical(dim(d), c(length(x), length(x)), info = paste(distance, "single-arg"))

          d2 <- proxy::dist(x, x, method = distance, window.size = 15L)

          expect_equal(d2, d, tolerance = 0, check.attributes = FALSE,
                       info = paste(distance, "double-arg"))

          d3 <- proxy::dist(x[1L], x, method = distance, window.size = 15L)
          class(d3) <- "matrix"

          expect_identical(dim(d3), c(1L, length(x)), info = paste(distance, "one-vs-many"))

          d4 <- proxy::dist(x, x[1L], method = distance, window.size = 15L)
          class(d4) <- "matrix"

          expect_identical(dim(d4), c(length(x), 1L), info = paste(distance, "many-vs-one"))

          ## dtw_lb will give different results below because of how it works
          if (distance == "dtw_lb") next

          expect_equal(d3, d[1L, , drop = FALSE], tolerance = 0, check.attributes = FALSE,
                       info = paste(distance, "one-vs-many-vs-distmat"))

          expect_equal(d4, d[ , 1L, drop = FALSE], tolerance = 0, check.attributes = FALSE,
                       info = paste(distance, "many-vs-one-vs-distmat"))
     }
})

# =================================================================================================
# proxy pairwise distances
# =================================================================================================

test_that("Included proxy distances can be called for pairwise = TRUE and give expected length", {
     for (distance in included_distances) {
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

          if (identical(Sys.getenv("NOT_CRAN"), "true")) {
               d3 <- proxy::dist(x[1L:3L], x[4L:6L], method = distance, window.size = 15L, pairwise = TRUE)
               expect_equal_to_reference(d3, paste0("rds/pdist_", distance, ".rds"))
          }
     }
})

# =================================================================================================
# distance functions
# =================================================================================================

ctrl <- new("dtwclustControl", window.size = 18L)
x <- data_reinterpolated[1L:20L]
centroids <- x[c(1L, 15L)]
attr(centroids, "id_cent") <- c(1L, 15L)

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

test_that("Errors in distance argument are correctly detected.", {
     expect_error(dtwclust(data_matrix, k = 20, distance = mean), "proxy", info = "Function")

     expect_error(dtwclust(data_matrix, k = 20, distance = NULL), "proxy", info = "NULL")

     expect_error(dtwclust(data_matrix, k = 20, distance = NA), "proxy", info = "NA")

     expect_error(dtwclust(data_matrix, k = 20, distance = "dummy"), "proxy", info = "Unregistered")

     expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                  "different length", info = "LBK")

     expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                  "different length", info = "LBI")

     expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                  "different length", info = "DTW_LB")
})
