context("Test pairwise")

# =================================================================================================
# sbd
# =================================================================================================

pdist_sbd <- proxy::dist(data[1:50], data[51:100], method = "SBD", force.pairwise = TRUE)

test_that("pairwise sbd distance gives the same result as reference",
          my_expect_equal_to_reference(pdist_sbd))

# =================================================================================================
# lbk
# =================================================================================================

pdist_lbk <- proxy::dist(data_list[1:50], data_list[51:100], method = "LBK",
                         window.size = 18L, force.pairwise = TRUE)

test_that("pairwise lbk distance gives the same result as reference",
          my_expect_equal_to_reference(pdist_lbk))

# =================================================================================================
# lbi
# =================================================================================================

pdist_lbi <- proxy::dist(data_list[1:50], data_list[51:100], method = "LBI",
                         window.size = 18L, force.pairwise = TRUE)

test_that("pairwise lbi distance gives the same result as reference",
          my_expect_equal_to_reference(pdist_lbi))

# =================================================================================================
# dtw_lb
# =================================================================================================

pdist_dtw_lb <- proxy::dist(data_list[1:50], data_list[51:100], method = "DTW_LB",
                            window.size = 18L, force.pairwise = TRUE)

test_that("pairwise dtw_lb distance gives the same result as reference",
          my_expect_equal_to_reference(pdist_dtw_lb))
