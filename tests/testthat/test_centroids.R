context("Test centroids")

ctrl <- new("dtwclustControl", window.size = 18L)
x <- data_reinterpolated[1L:20L]
k <- 2L
cl_id <- rep(c(1L, 2L), each = length(x) / 2L)

# =================================================================================================
# mean
# =================================================================================================

test_that("Operations with mean centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "mean")

     pc_mean <- family@allcent(x,
                               cl_id = cl_id,
                               k = k,
                               cent = x[c(1L,20L)],
                               cl_old = 0L)

     expect_identical(length(pc_mean), k)

     ## ---------------------------------------------------------- multivariate
     mv_mean <- family@allcent(data_multivariate,
                               cl_id = cl_id,
                               k = k,
                               cent = data_multivariate[c(1L,20L)],
                               cl_old = 0L)

     expect_identical(length(mv_mean), k)

     expect_identical(dim(mv_mean[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(pc_mean, file_name(pc_mean), info = "Univariate")
     expect_equal_to_reference(mv_mean, file_name(mv_mean), info = "Multivariate")
})

# =================================================================================================
# median
# =================================================================================================

test_that("Operations with median centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "median")

     pc_median <- family@allcent(x,
                                 cl_id = cl_id,
                                 k = k,
                                 cent = x[c(1L,20L)],
                                 cl_old = 0L)

     expect_identical(length(pc_median), k)

     ## ---------------------------------------------------------- multivariate
     mv_median <- family@allcent(data_multivariate,
                                 cl_id = cl_id,
                                 k = k,
                                 cent = data_multivariate[c(1L,20L)],
                                 cl_old = 0L)

     expect_identical(length(mv_median), k)

     expect_identical(dim(mv_median[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(pc_median, file_name(pc_median), info = "Univariate")
     expect_equal_to_reference(mv_median, file_name(mv_median), info = "Multivariate")
})

# =================================================================================================
# shape
# =================================================================================================

test_that("Operations with shape centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "shape")

     pc_shape <- family@allcent(x,
                                cl_id = cl_id,
                                k = k,
                                cent = x[c(1L,20L)],
                                cl_old = 0L)

     expect_identical(length(pc_shape), k)

     ## ---------------------------------------------------------- multivariate
     mv_shape <- family@allcent(data_multivariate,
                                cl_id = cl_id,
                                k = k,
                                cent = data_multivariate[c(1L,20L)],
                                cl_old = 0L)

     expect_identical(length(mv_shape), k)

     expect_identical(dim(mv_shape[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(pc_shape, file_name(pc_shape), info = "Univariate")
     expect_equal_to_reference(mv_shape, file_name(mv_shape), info = "Multivariate")
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

     pc_pam <- family@allcent(x,
                              cl_id = cl_id,
                              k = k,
                              cent = x[c(1L,20L)],
                              cl_old = 0L)

     expect_identical(length(pc_pam), k)

     expect_null(as.list(environment(family@allcent))$distmat)

     ## ---------------------------------------------------------- univariate with distmat
     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "sbd",
                   allcent = "pam",
                   distmat = proxy::dist(x, method = "sbd"))

     pc_pam_distmat <- family@allcent(x,
                                      cl_id = cl_id,
                                      k = k,
                                      cent = x[c(1L,20L)],
                                      cl_old = 0L)

     expect_identical(length(pc_pam_distmat), k)

     expect_false(is.null(as.list(environment(family@allcent))$distmat))

     expect_identical(pc_pam, pc_pam_distmat)

     ## ---------------------------------------------------------- multivariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   dist = "dtw_basic",
                   allcent = "pam")

     mv_pam <- family@allcent(data_multivariate,
                              cl_id = cl_id,
                              k = k,
                              cent = data_multivariate[c(1L,20L)],
                              cl_old = 0L)

     expect_identical(length(mv_pam), k)

     expect_identical(dim(mv_pam[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(pc_pam, file_name(pc_pam), info = "Univariate without distmat")
     expect_equal_to_reference(mv_pam, file_name(mv_pam), info = "Multivariate")
})

# =================================================================================================
# dba
# =================================================================================================

test_that("Operations with dba centroid give same results as references.", {
     ## ---------------------------------------------------------- univariate
     family <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "dba")

     pc_dba <- family@allcent(x,
                              cl_id = cl_id,
                              k = k,
                              cent = x[c(1L,20L)],
                              cl_old = 0L)

     expect_identical(length(pc_dba), k)

     ## ---------------------------------------------------------- multivariate
     mv_dba <- family@allcent(data_multivariate,
                              cl_id = cl_id,
                              k = k,
                              cent = data_multivariate[c(1L,20L)],
                              cl_old = 0L)

     expect_identical(length(mv_dba), k)

     expect_identical(dim(mv_dba[[1L]]), dim(data_multivariate[[1L]]))

     skip_on_cran()

     expect_equal_to_reference(pc_dba, file_name(pc_dba), info = "Univariate")
     expect_equal_to_reference(mv_dba, file_name(mv_dba), info = "Multivariate")
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

     pc_colMeans <- dtwclust(data_matrix, type = "partitional", k = 20,
                             distance = "sbd", centroid = mycent,
                             preproc = NULL, control = ctrl, seed = 123)

     pc_colMeans <- reset_nondeterministic(pc_colMeans)

     skip_on_cran()

     expect_equal_to_reference(pc_colMeans, file_name(pc_colMeans), info = "Custom colMeans")
})

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
})
