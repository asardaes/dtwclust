context("Test fuzzy")

# Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat) {
     lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# =================================================================================================
# Valid input
# =================================================================================================

fc <- dtwclust(data_subset, type = "fuzzy", k = 4:6,
               preproc = acf_fun, distance = "L2",
               seed = 123)

fc <- lapply(fc, reset_nondeterministic)

test_that("Fuzzy clustering gives the same result as reference",
          my_expect_equal_to_reference(fc))

# =================================================================================================
# Invalid input
# =================================================================================================

test_that("Fuzzy clustering doesn't allow series with different lengths",
          expect_error(dtwclust(data_subset, type = "fuzzy", k = 4),
                       "different length"))
