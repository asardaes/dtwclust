context("Test repetitions")

ctrl@nrep <- 4L

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

ctrl@pam.precompute <- FALSE

rep_npampre <- dtwclust(data_matrix, type = "partitional", k = 20:21,
                       distance = "L1", centroid = "pam",
                       preproc = NULL, control = ctrl, seed = 321)

rep_npampre <- lapply(rep_npampre, function(x) {
     x <- reset_nondeterministic(x)
     x@control@pam.precompute <- TRUE
     x
})

rep_pampre <- lapply(rep_pampre, function(x) {
     x@distmat <- NULL
     x
})

test_that("Repetitions without pam.precompute give appropriate length",
          expect_identical(length(rep_npampre), ctrl@nrep * 2L))

test_that("Repetitions without pam.precompute give the same result as with pam.precompute (adjusted distmat)",
          expect_identical(rep_pampre, rep_npampre))

ctrl@nrep <- 1L
ctrl@pam.precompute <- TRUE
