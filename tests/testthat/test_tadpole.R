context("Test TADPole")

# =================================================================================================
# Valid input
# =================================================================================================

pc_tadp <- dtwclust(data_list[1:50], type = "tadpole", k = 10, dc = 1.5, control = ctrl)

pc_tadp <- reset_nondeterministic(pc_tadp)

test_that("TADPole clustering gives the same result as reference",
          my_expect_equal_to_reference(pc_tadp))

# =================================================================================================
# Invalid input
# =================================================================================================

test_that("TADPole clustering doesn't allow series with different lengths",
          expect_error(dtwclust(data, type = "tadpole", k = 20, dc = 1.5, control = ctrl),
                       "same length"))
