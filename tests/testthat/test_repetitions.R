context("Test repetitions")

ctrl@nrep <- 8L

# =================================================================================================
# Precomputing distmat
# =================================================================================================

rep_pampre <- dtwclust(data_matrix, type = "partitional", k = 20:21,
                       distance = "L1", centroid = "pam",
                       preproc = NULL, control = ctrl, seed = 321)

rep_pampre <- lapply(rep_pampre, reset_nondeterministic)

test_that("Repetitions with pam.precompute give appropriate length",
          expect_identical(length(rep_pampre), ctrl@nrep * 2L))

test_that("Repetitions with pam.precompute give the same result as reference",
          my_expect_equal_to_reference(rep_pampre))

# =================================================================================================
# Without precomputing distmat
# =================================================================================================

rep_median <- dtwclust(data_matrix, type = "partitional", k = 19:20,
                       distance = "L1", centroid = "median",
                       preproc = NULL, control = ctrl, seed = 321)

rep_median <- lapply(rep_median, reset_nondeterministic)

test_that("Repetitions without pam.precompute give appropriate length",
          expect_identical(length(rep_median), ctrl@nrep * 2L))

test_that("Repetitions without pam.precompute give the same result as reference",
          my_expect_equal_to_reference(rep_median))

ctrl@nrep <- 1L
