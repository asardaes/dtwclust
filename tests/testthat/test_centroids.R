context("Test centroids")

# =================================================================================================
# mean
# =================================================================================================

## Will not converge
suppressWarnings(
     pc_mean <- dtwclust(data_matrix, type = "partitional", k = 20,
                         distance = "sbd", centroid = "mean",
                         preproc = NULL, control = ctrl, seed = 123)
)

pc_mean <- reset_nondeterministic(pc_mean)

test_that("mean centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_mean))

# =================================================================================================
# median
# =================================================================================================

pc_median <- dtwclust(data_matrix, type = "partitional", k = 20,
                      distance = "sbd", centroid = "median",
                      preproc = NULL, control = ctrl, seed = 123)

pc_median <- reset_nondeterministic(pc_median)

test_that("median centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_median))

# =================================================================================================
# shape
# =================================================================================================

pc_shape <- dtwclust(data, type = "partitional", k = 20,
                     distance = "sbd", centroid = "shape",
                     preproc = NULL, control = ctrl, seed = 123)

pc_shape <- reset_nondeterministic(pc_shape)

test_that("shape centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_shape))

# =================================================================================================
# pam
# =================================================================================================

pc_pam <- dtwclust(data, type = "partitional", k = 20,
                   distance = "sbd", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

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
                        preproc = NULL, control = ctrl, seed = 123)
)

pc_dba <- reset_nondeterministic(pc_dba)

test_that("dba centroid gives the same result as reference",
          my_expect_equal_to_reference(pc_dba))

# =================================================================================================
# colMeans
# =================================================================================================

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
