context("Test hierarchical")

# =================================================================================================
# all
# =================================================================================================

hc_all <- dtwclust(data, type = "hierarchical", k = 20:21,
                   distance = "sbd", method = "all",
                   control = ctrl)

hc_all <- lapply(hc_all, reset_nondeterministic)

test_that("hierarchical clustering gives the same result as reference",
          my_expect_equal_to_reference(hc_all, TRUE))

# =================================================================================================
# non-symmetric distance
# =================================================================================================

hc_lbi <- dtwclust(data_list, type = "hierarchical", k = 19:20,
                   distance = "lbi", method = "all",
                   control = ctrl)

hc_lbi <- lapply(hc_lbi, reset_nondeterministic)

test_that("hierarchical clustering with non-symmetric distance gives the same result as reference",
          my_expect_equal_to_reference(hc_lbi))

# =================================================================================================
# distance function
# =================================================================================================

test_that("distance function gives error",
          expect_error(dtwclust(data_matrix, k = 20, type = "hierarchical", distance = mean),
                       "proxy"))

# =================================================================================================
# unregistered distance
# =================================================================================================

test_that("unregistered distance gives error",
          expect_error(dtwclust(data_matrix, k = 20, type = "hierarchical", distance = "dummy"),
                       "proxy"))
