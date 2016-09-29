context("Test fuzzy")

# Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat, ...) {
     lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# =================================================================================================
# invalid combinations
# =================================================================================================

test_that("Invalid combinations in hierarchical clustering are detected.", {
     expect_error(dtwclust(data_matrix, type = "f", k = 101L), "more clusters")

     ctrl <- list(window.size = 15L)

     expect_error(dtwclust(data, type = "f", distance = "lbk", control = ctrl), "different length")
     expect_error(dtwclust(data, type = "f", distance = "lbi", control = ctrl), "different length")
     expect_error(dtwclust(data, type = "f", distance = "dtw_lb", control = ctrl), "different length")
     expect_error(dtwclust(data, type = "f", distance = "dtw"), "different length")
     expect_error(dtwclust(data, type = "f", distance = "dtw2"), "different length")
     expect_error(dtwclust(data, type = "f", distance = "dtw_basic"), "different length")
     expect_error(dtwclust(data, type = "f", distance = "sbd"), "different length")

     expect_error(dtwclust(data, type = "f", preproc = "zscore"), "preprocessing")
     expect_error(dtwclust(data, type = "f", preproc = reinterpolate), "preprocessing.*arguments")

     expect_error(dtwclust(data_matrix, type = "f", distance = mean), "proxy", info = "Function")
     expect_error(dtwclust(data_matrix, type = "f", distance = NULL), "proxy", info = "NULL")
     expect_error(dtwclust(data_matrix, type = "f", distance = NA), "proxy", info = "NA")
     expect_error(dtwclust(data_matrix, type = "f", distance = "dummy"), "proxy", info = "Unregistered")

     expect_warning(dtwclust(data_subset, type = "f",
                             preproc = acf_fun, distance = "L2",
                             centroid = "mean"),
                    "centroid.*fcm")

     expect_warning(dtwclust(data_subset, type = "f",
                             preproc = acf_fun, distance = "L2",
                             centroid = "median"),
                    "centroid.*fcm")

     expect_warning(dtwclust(data_subset, type = "f",
                             preproc = acf_fun, distance = "L2",
                             centroid = "shape"),
                    "centroid.*fcm")

     expect_warning(dtwclust(data_subset, type = "f",
                             preproc = acf_fun, distance = "L2",
                             centroid = "dba"),
                    "centroid.*fcm")

     expect_warning(dtwclust(data_subset, type = "f",
                             preproc = acf_fun, distance = "L2",
                             centroid = "pam"),
                    "centroid.*fcm")
})

# =================================================================================================
# multiple k
# =================================================================================================

test_that("Multiple k works as expected.", {
     fc_k <- dtwclust(data_subset, type = "f", k = 2L:5L,
                      preproc = acf_fun, distance = "L2",
                      seed = 123)

     expect_identical(length(fc_k), 4L)

     skip_on_cran()

     fc_k <- lapply(fc_k, reset_nondeterministic)
     expect_equal_to_reference(fc_k, file_name(fc_k))
})

# =================================================================================================
# valid input
# =================================================================================================

test_that("Fuzzy clustering works as expected.", {
     skip_on_cran()

     ## ---------------------------------------------------------- univariate
     fc <- dtwclust(data_subset, type = "fuzzy", k = 4L,
                    preproc = acf_fun, distance = "L2",
                    seed = 123)

     fc <- reset_nondeterministic(fc)

     expect_equal_to_reference(fc, file_name(fc))

     ## ---------------------------------------------------------- multivariate
     fc_mv <- dtwclust(data_multivariate[1L:20L], type = "fuzzy", k = 4L,
                    preproc = acf_fun, distance = "dtw_basic",
                    seed = 123)

     fc_mv <- reset_nondeterministic(fc_mv)

     expect_equal_to_reference(fc_mv, file_name(fc_mv))
})
