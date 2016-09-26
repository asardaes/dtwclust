context("Test hierarchical")

# =================================================================================================
# invalid combinations
# =================================================================================================

test_that("Invalid combinations in hierarchical clustering are detected.", {
     expect_error(dtwclust(data, type = "h", k = 101L), "more clusters")

     expect_error(dtwclust(data, type = "h", distance = "lbk"), "different length")
     expect_error(dtwclust(data, type = "h", distance = "lbi"), "different length")
     expect_error(dtwclust(data, type = "h", distance = "dtw_lb"), "different length")

     expect_error(dtwclust(data, type = "h", preproc = "zscore"), "preprocessing")

     expect_error(dtwclust(data_matrix, type = "h", distance = mean), "proxy", info = "Function")
     expect_error(dtwclust(data_matrix, type = "h", distance = NULL), "proxy", info = "NULL")
     expect_error(dtwclust(data_matrix, type = "h", distance = NA), "proxy", info = "NA")
     expect_error(dtwclust(data_matrix, type = "h", distance = "dummy"), "proxy", info = "Unregistered")
})

# =================================================================================================
# multiple k
# =================================================================================================

test_that("Multiple k works as expected.", {
     hc_k <- dtwclust(data_reinterpolated, type = "h", k = 2L:5L,
                      distance = "L2", seed = 938)

     expect_identical(length(hc_k), 4L)

     skip_on_cran()

     hc_k <- lapply(hc_k, reset_nondeterministic)
     expect_equal_to_reference(hc_k, file_name(hc_k))
})

# =================================================================================================
# hierarchical algorithms
# =================================================================================================

test_that("Hierarchical clustering works as expected.", {
     skip_on_cran()

     ## ---------------------------------------------------------- all
     hc_all <- dtwclust(data, type = "hierarchical", k = 20L,
                        distance = "sbd", method = "all")

     hc_all <- lapply(hc_all, reset_nondeterministic)

     expect_equal_to_reference(hc_all, file_name(hc_all))

     ## ---------------------------------------------------------- non-symmetric
     hc_lbi <- dtwclust(data_reinterpolated, type = "hierarchical", k = 20L,
                        distance = "lbi", method = "all",
                        control = list(window.size = 17L))

     hc_lbi <- lapply(hc_lbi, reset_nondeterministic)

     expect_equal_to_reference(hc_lbi, file_name(hc_lbi))

     ## ---------------------------------------------------------- custom centroid
     hc_cent <- dtwclust(data, type = "hierarchical", k = 20L,
                         distance = "sbd", method = "all",
                         preproc = zscore, centroid = shape_extraction,
                         seed = 320)

     hc_cent <- lapply(hc_cent, reset_nondeterministic)

     expect_equal_to_reference(hc_cent, file_name(hc_cent))
})
