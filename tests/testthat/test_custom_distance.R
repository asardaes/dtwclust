context("Test custom distance")

# =================================================================================================
# Registered with proxy
# =================================================================================================

ndtw <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = symmetric2,
              distance.only = TRUE, ...)$normalizedDistance
}

if (!pr_DB$entry_exists("nDTW"))
     proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "Normalized DTW with L1 norm")

pc_ndtw <- dtwclust(data_list[1:20], k = 4, distance = "nDTW",
                    control = ctrl, seed = 8319)

pc_ndtw <- reset_nondeterministic(pc_ndtw)

test_that("custom distance gives the same result as reference",
          my_expect_equal_to_reference(pc_ndtw))

# =================================================================================================
# A symmetric computation of the above (since lengths are equal)
# =================================================================================================

ctrl@symmetric <- TRUE

pc_ndtw <- dtwclust(data_list[1:20], k = 4, distance = "nDTW",
                    control = ctrl, seed = 8319)

pc_ndtw <- reset_nondeterministic(pc_ndtw)
pc_ndtw@control@symmetric <- FALSE # just for the check below

test_that("custom SYMMETRIC distance gives the same result as NON-SYMMETRIC when appropriate",
          my_expect_equal_to_reference(pc_ndtw))

ctrl@symmetric <- FALSE

# =================================================================================================
# Registered with proxy and with custom arguments
# =================================================================================================

ndtw2 <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = asymmetric,
              distance.only = TRUE, ...)$normalizedDistance
}

if (!pr_DB$entry_exists("nDTW2"))
     proxy::pr_DB$set_entry(FUN = ndtw2, names=c("nDTW2"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "Normalized asymmetric DTW with L1 norm")

pc_ndtw2 <- dtwclust(data_subset, k = 4, distance = "nDTW2",
                     control = ctrl, seed = 8319,
                     open.begin = TRUE, open.end = TRUE)

pc_ndtw2 <- reset_nondeterministic(pc_ndtw2)

test_that("custom distance with custom arguments gives the same result as reference",
          my_expect_equal_to_reference(pc_ndtw2))
