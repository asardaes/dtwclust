context("Test centroids and CVIs")

## CVIs require data and must be evaluated before resetting nondeterministic values,
## since the functions in the family slot are needed
ctrl@save.data <- TRUE

# =================================================================================================
# mean
# =================================================================================================

## Will not converge
suppressWarnings(
     pc_mean <- dtwclust(data_matrix, type = "partitional", k = 20,
                         distance = "sbd", centroid = "mean",
                         preproc = NULL, control = ctrl, seed = 123)
)

cvi_mean <- cvi(pc_mean, labels)

test_that("mean centroid CVI gives the same result as reference",
          my_expect_equal_to_reference(cvi_mean))

pc_mean <- reset_nondeterministic(pc_mean)

test_that("mean centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_mean))

# =================================================================================================
# median
# =================================================================================================

pc_median <- dtwclust(data_matrix, type = "partitional", k = 20,
                      distance = "sbd", centroid = "median",
                      preproc = NULL, control = ctrl, seed = 123)

cvi_median <- cvi(pc_median, labels)

test_that("median centroid CVI gives the same result as reference",
          my_expect_equal_to_reference(cvi_median))

pc_median <- reset_nondeterministic(pc_median)

test_that("median centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_median))

# =================================================================================================
# shape
# =================================================================================================

pc_shape <- dtwclust(data, type = "partitional", k = 20,
                     distance = "sbd", centroid = "shape",
                     preproc = NULL, control = ctrl, seed = 123)

cvi_shape <- cvi(pc_shape, labels)

test_that("shape centroid CVI gives the same result as reference",
          my_expect_equal_to_reference(cvi_shape))

pc_shape <- reset_nondeterministic(pc_shape)

test_that("shape centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_shape))

# =================================================================================================
# pam
# =================================================================================================

pc_pam <- dtwclust(data, type = "partitional", k = 20,
                   distance = "sbd", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

cvi_pam <- cvi(pc_pam, labels)

test_that("pam centroid CVI gives the same result as reference",
          my_expect_equal_to_reference(cvi_pam))

pc_pam <- reset_nondeterministic(pc_pam)

test_that("pam centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_pam, TRUE))

# =================================================================================================
# dba
# =================================================================================================

## Will not converge
suppressWarnings(
     pc_dba <- dtwclust(data_subset, type = "partitional", k = 4,
                        distance = "sbd", centroid = "dba",
                        preproc = NULL, control = ctrl, seed = 123,
                        dba.alignment = "dtw")
)

suppressWarnings(cvi_dba <- cvi(pc_dba, labels_subset))

test_that("dba centroid CVI gives the same result as reference",
          my_expect_equal_to_reference(cvi_dba))

pc_dba <- reset_nondeterministic(pc_dba)

test_that("dba centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_dba))

suppressWarnings(
     pc_dba_basic <- dtwclust(data_subset, type = "partitional", k = 4,
                              distance = "sbd", centroid = "dba",
                              preproc = NULL, control = ctrl, seed = 123,
                              dba.alignment = "dtw_basic")
)

pc_dba_basic <- reset_nondeterministic(pc_dba_basic)

test_that("dba_basic centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_dba_basic))

# =================================================================================================
# colMeans
# =================================================================================================

ctrl@save.data <- FALSE

mycent <- function(x, cl_id, k, cent, cl_old, ...) {
     x_split <- split(x, cl_id)

     x_split <- lapply(x_split, function(xx) do.call(rbind, xx))

     new_cent <- lapply(x_split, colMeans)

     new_cent
}

## Will not converge
suppressWarnings(
     pc_colMeans <- dtwclust(data_matrix, type = "partitional", k = 20,
                             distance = "sbd", centroid = mycent,
                             preproc = NULL, control = ctrl, seed = 123)
)

pc_colMeans <- reset_nondeterministic(pc_colMeans)

test_that("custom centroid function gives the same result as reference",
          my_expect_equal_to_reference(pc_colMeans))

# =================================================================================================
# invalid centroid
# =================================================================================================

test_that("invalid centroid for series with different lengths gives error",
          expect_error(dtwclust(data, k = 20, distance = "sbd", centroid = "mean"),
                       "different length"))
# =================================================================================================
# mean multivariate
# =================================================================================================

mv_mean <- dtwclust(data_multivariate, type = "partitional", k = 4,
                    distance = "dtw", centroid = "mean",
                    preproc = NULL, control = ctrl, seed = 123,
                    dist.method = "L1")

mv_mean <- reset_nondeterministic(mv_mean)

test_that("multivariate mean centroid gives the same result as reference",
          my_expect_equal_to_reference(mv_mean))

# =================================================================================================
# median multivariate
# =================================================================================================

mv_median <- dtwclust(data_multivariate, type = "partitional", k = 4,
                      distance = "dtw2", centroid = "median",
                      preproc = NULL, control = ctrl, seed = 123)

mv_median <- reset_nondeterministic(mv_median)

test_that("multivariate median centroid gives the same result as reference",
          my_expect_equal_to_reference(mv_median))

# =================================================================================================
# shape multivariate
# =================================================================================================

mv_shape <- dtwclust(data_multivariate, type = "partitional", k = 4,
                     distance = "dtw", centroid = "shape",
                     preproc = NULL, control = ctrl, seed = 123,
                     dist.method = "L1")

mv_shape <- reset_nondeterministic(mv_shape)

test_that("multivariate shape centroid gives the same result as reference",
          my_expect_equal_to_reference(mv_shape))

# =================================================================================================
# pam multivariate
# =================================================================================================

mv_pam <- dtwclust(data_multivariate, type = "partitional", k = 4,
                   distance = "dtw2", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

mv_pam <- reset_nondeterministic(mv_pam)

test_that("multivariate pam centroid gives the same result as reference",
          my_expect_equal_to_reference(mv_pam, TRUE))

# =================================================================================================
# dba multivariate
# =================================================================================================

## Will not converge
suppressWarnings(
     mv_dba <- dtwclust(data_multivariate, type = "partitional", k = 4,
                        distance = "dtw", centroid = "dba",
                        preproc = NULL, control = ctrl, seed = 123,
                        dist.method = "L1")
)

mv_dba <- reset_nondeterministic(mv_dba)

test_that("multivariate dba centroid gives the same result as reference",
          my_expect_equal_to_reference(mv_dba))
