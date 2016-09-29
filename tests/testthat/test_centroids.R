context("Test included centroids")

ctrl <- new("dtwclustControl", window.size = 18L)
x <- data_reinterpolated[1L:20L]
k <- 2L
cl_id <- rep(c(1L, 2L), each = length(x) / 2L)

# =================================================================================================
# invalid centroid
# =================================================================================================

test_that("Errors in centroid argument are correctly detected.", {
     expect_error(dtwclust(data, centroid = "mean"),
                  "different length", "mean")

     expect_error(dtwclust(data, centroid = "mean"),
                  "different length", "median")

     expect_error(dtwclust(data, centroid = "fcm"),
                  "arg", info = "FCM is for fuzzy")

     expect_error(dtwclust(data, centroid = "goku"),
                  "arg", info = "Unknown centroid")

     expect_error(dtwclust(data, centroid = NULL),
                  "definition", info = "NULL centroid")

     expect_error(dtwclust(data, centroid = NA),
                  "definition", info = "NA centroid")

     expect_error(dtwclust(data, centroid = function(x, ...) { x }),
                  "centroid.*arguments")
})

# =================================================================================================
# mean
# =================================================================================================

test_that("Operations with mean centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "mean")

     cent_mean <- family@allcent(x,
                                 cl_id = cl_id,
                                 k = k,
                                 cent = x[c(1L,20L)],
                                 cl_old = 0L)

     expect_identical(length(cent_mean), k)

     ## ---------------------------------------------------------- multivariate
     cent_mv_mean <- family@allcent(data_multivariate,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = data_multivariate[c(1L,20L)],
                                    cl_old = 0L)

     expect_identical(length(cent_mv_mean), k)

     expect_identical(dim(cent_mv_mean[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(cent_mean, file_name(cent_mean), info = "Univariate")
     expect_equal_to_reference(cent_mv_mean, file_name(cent_mv_mean), info = "Multivariate")
})

# =================================================================================================
# median
# =================================================================================================

test_that("Operations with median centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "median")

     cent_median <- family@allcent(x,
                                   cl_id = cl_id,
                                   k = k,
                                   cent = x[c(1L,20L)],
                                   cl_old = 0L)

     expect_identical(length(cent_median), k)

     ## ---------------------------------------------------------- multivariate
     cent_mv_median <- family@allcent(data_multivariate,
                                      cl_id = cl_id,
                                      k = k,
                                      cent = data_multivariate[c(1L,20L)],
                                      cl_old = 0L)

     expect_identical(length(cent_mv_median), k)

     expect_identical(dim(cent_mv_median[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(cent_median, file_name(cent_median), info = "Univariate")
     expect_equal_to_reference(cent_mv_median, file_name(cent_mv_median), info = "Multivariate")
})

# =================================================================================================
# shape
# =================================================================================================

test_that("Operations with shape centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "shape")

     cent_shape <- family@allcent(x,
                                  cl_id = cl_id,
                                  k = k,
                                  cent = x[c(1L,20L)],
                                  cl_old = 0L)

     expect_identical(length(cent_shape), k)

     ## ---------------------------------------------------------- multivariate
     cent_mv_shape <- family@allcent(data_multivariate,
                                     cl_id = cl_id,
                                     k = k,
                                     cent = data_multivariate[c(1L,20L)],
                                     cl_old = 0L)

     expect_identical(length(cent_mv_shape), k)

     expect_identical(dim(cent_mv_shape[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(cent_shape, file_name(cent_shape), info = "Univariate")
     expect_equal_to_reference(cent_mv_shape, file_name(cent_mv_shape), info = "Multivariate")
})

# =================================================================================================
# pam
# =================================================================================================

test_that("Operations with pam centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate without distmat
     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "sbd",
                   allcent = "pam")

     cent_pam <- family@allcent(x,
                                cl_id = cl_id,
                                k = k,
                                cent = x[c(1L,20L)],
                                cl_old = 0L)

     expect_identical(length(cent_pam), k)

     expect_null(as.list(environment(family@allcent))$distmat)

     ## ---------------------------------------------------------- univariate with distmat
     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "sbd",
                   allcent = "pam",
                   distmat = proxy::dist(x, method = "sbd"))

     cent_pam_distmat <- family@allcent(x,
                                        cl_id = cl_id,
                                        k = k,
                                        cent = x[c(1L,20L)],
                                        cl_old = 0L)

     expect_identical(length(cent_pam_distmat), k)

     expect_false(is.null(as.list(environment(family@allcent))$distmat))

     expect_identical(cent_pam, cent_pam_distmat)

     ## ---------------------------------------------------------- multivariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "dtw_basic",
                   allcent = "pam")

     cent_mv_pam <- family@allcent(data_multivariate,
                                   cl_id = cl_id,
                                   k = k,
                                   cent = data_multivariate[c(1L,20L)],
                                   cl_old = 0L)

     expect_identical(length(cent_mv_pam), k)

     expect_identical(dim(cent_mv_pam[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(cent_pam, file_name(cent_pam), info = "Univariate without distmat")
     expect_equal_to_reference(cent_mv_pam, file_name(cent_mv_pam), info = "Multivariate")
})

# =================================================================================================
# dba
# =================================================================================================

test_that("Operations with dba centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "dba")

     cent_dba <- family@allcent(x,
                                cl_id = cl_id,
                                k = k,
                                cent = x[c(1L,20L)],
                                cl_old = 0L)

     expect_identical(length(cent_dba), k)

     ## ---------------------------------------------------------- multivariate
     cent_mv_dba <- family@allcent(data_multivariate,
                                   cl_id = cl_id,
                                   k = k,
                                   cent = data_multivariate[c(1L,20L)],
                                   cl_old = 0L)

     expect_identical(length(cent_mv_dba), k)

     expect_identical(dim(cent_mv_dba[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(cent_dba, file_name(cent_dba), info = "Univariate")
     expect_equal_to_reference(cent_mv_dba, file_name(cent_mv_dba), info = "Multivariate")
})

# =================================================================================================
# custom
# =================================================================================================

test_that("Operations with custom centroid give same results as references.", {
     mycent <- function(x, cl_id, k, cent, cl_old, ...) {
          x_split <- split(x, cl_id)

          x_split <- lapply(x_split, function(xx) do.call(rbind, xx))

          new_cent <- lapply(x_split, colMeans)

          new_cent
     }

     cent_colMeans <- dtwclust(data_matrix, type = "partitional", k = 20,
                               distance = "sbd", centroid = mycent,
                               preproc = NULL, control = ctrl, seed = 123)

     cent_colMeans <- reset_nondeterministic(cent_colMeans)

     skip_on_cran()

     expect_equal_to_reference(cent_colMeans, file_name(cent_colMeans), info = "Custom colMeans")
})
