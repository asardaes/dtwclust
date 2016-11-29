context("\tFuzzy")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

## Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat, ...) {
    lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# =================================================================================================
# multiple k
# =================================================================================================

test_that("Multiple k works as expected.", {
    fc_k <- dtwclust(data_subset, type = "f", k = 2L:5L,
                     preproc = acf_fun, distance = "L2",
                     seed = 123)

    expect_identical(length(fc_k), 4L)

    fc_k <- lapply(fc_k, reset_nondeterministic)
    assign("fc_k", fc_k, persistent)
})

# =================================================================================================
# valid input
# =================================================================================================

test_that("Fuzzy clustering works as expected.", {
    ## ---------------------------------------------------------- univariate
    fc <- dtwclust(data_subset, type = "fuzzy", k = 4L,
                   preproc = acf_fun, distance = "L2",
                   seed = 123)

    fc <- reset_nondeterministic(fc)

    assign("fc", fc, persistent)

    ## ---------------------------------------------------------- multivariate
    dmv <- reinterpolate(data_multivariate, new.length = max(sapply(data_multivariate, NROW)))

    fc_mv <- dtwclust(dmv, type = "fuzzy", k = 4L,
                      distance = "dtw_basic",
                      seed = 123)

    fc_mv <- reset_nondeterministic(fc_mv)

    assign("fc_mv", fc_mv, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
