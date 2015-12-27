context("Test hierarchical")

# =================================================================================================
# all
# =================================================================================================

hc_all <- dtwclust(data, type = "hierarchical", k = 20,
                   distance = "sbd", method = "all")

hc_all <- lapply(hc_all, reset_nondeterministic)

test_that("hierarchical clustering gives the same result as reference",
          expect_equal_to_reference(hc_all, "hc_all.rds"))

# =================================================================================================
# non-symmetric distance
# =================================================================================================

hc_lbi <- dtwclust(data_list, type = "hierarchical", k = 20,
                   distance = "lbi", method = "all",
                   control = list(window.size = 18L))

hc_lbi <- lapply(hc_lbi, reset_nondeterministic)

test_that("hierarchical clustering with non-symmetric distance gives the same result as reference",
          expect_equal_to_reference(hc_lbi, "hc_lbi.rds"))

# =================================================================================================
# distance function
# =================================================================================================

test_that("distance function gives error",
          expect_error(dtwclust(data_matrix, k = 20, type = "hierarchical", distance = mean)))

# =================================================================================================
# unregistered distance
# =================================================================================================

test_that("unregistered distance gives error",
          expect_error(dtwclust(data_matrix, k = 20, type = "hierarchical", distance = "dummy")))
