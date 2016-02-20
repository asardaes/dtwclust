context("Test distances")

# =================================================================================================
# L2
# =================================================================================================

pc_l2 <- dtwclust(data_matrix, type = "partitional", k = 20,
                  distance = "L2", centroid = "pam",
                  preproc = NULL, control = ctrl, seed = 123)

pc_l2 <- reset_nondeterministic(pc_l2)

test_that("L2 distance gives the same result as reference",
          my_expect_equal_to_reference(pc_l2))

# =================================================================================================
# lbk
# =================================================================================================

pc_lbk <- dtwclust(data_matrix, type = "partitional", k = 20,
                   distance = "lbk", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

pc_lbk <- reset_nondeterministic(pc_lbk)

test_that("lbk distance gives the same result as reference",
          my_expect_equal_to_reference(pc_lbk))

# =================================================================================================
# lbi
# =================================================================================================

pc_lbi <- dtwclust(data_matrix, type = "partitional", k = 20,
                   distance = "lbi", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

pc_lbi <- reset_nondeterministic(pc_lbi)

test_that("lbi distance gives the same result as reference",
          my_expect_equal_to_reference(pc_lbi))

# =================================================================================================
# sbd
# =================================================================================================

pc_sbd <- dtwclust(data_matrix, type = "partitional", k = 20,
                   distance = "sbd", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

pc_sbd <- reset_nondeterministic(pc_sbd)

test_that("sbd distance gives the same result as reference",
          my_expect_equal_to_reference(pc_sbd, TRUE))

# =================================================================================================
# dtw_lb
# =================================================================================================

ctrl@pam.precompute <- FALSE

pc_dtw_lb <- dtwclust(data_matrix[1:20, ], type = "partitional", k = 4,
                      distance = "dtw_lb", centroid = "pam",
                      preproc = NULL, control = ctrl, seed = 123)

pc_dtw_lb <- reset_nondeterministic(pc_dtw_lb)

test_that("dtw_lb distance gives the same result as reference",
          my_expect_equal_to_reference(pc_dtw_lb))

ctrl@pam.precompute <- TRUE

# =================================================================================================
# dtw
# =================================================================================================

pc_dtw <- dtwclust(data_subset, type = "partitional", k = 4,
                   distance = "dtw", centroid = "pam",
                   preproc = NULL, control = ctrl, seed = 123)

pc_dtw <- reset_nondeterministic(pc_dtw)

test_that("dtw distance gives the same result as reference",
          my_expect_equal_to_reference(pc_dtw))

# =================================================================================================
# dtw2
# =================================================================================================

pc_dtw2 <- dtwclust(data_subset, type = "partitional", k = 4,
                    distance = "dtw2", centroid = "pam",
                    preproc = NULL, control = ctrl, seed = 123)

pc_dtw2 <- reset_nondeterministic(pc_dtw2)

test_that("dtw2 distance gives the same result as reference",
          my_expect_equal_to_reference(pc_dtw2))

# =================================================================================================
# distance function
# =================================================================================================

test_that("distance function gives error",
          expect_error(dtwclust(data_matrix, k = 20, distance = mean),
                       "proxy"))

# =================================================================================================
# unregistered distance
# =================================================================================================

test_that("unregistered distance gives error",
          expect_error(dtwclust(data_matrix, k = 20, distance = "dummy"),
                       "proxy"))

# =================================================================================================
# invalid distance
# =================================================================================================

test_that("invalid distance for series with different lengths gives error",
          expect_error(dtwclust(data, k = 20, distance = "lbi",
                                control = list(window.size = 18L)),
                       "different length"))
