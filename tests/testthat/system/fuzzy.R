context("    Fuzzy")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

## Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat, ...) {
    lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# ==================================================================================================
# multiple k
# ==================================================================================================

test_that("Multiple k works as expected.", {
    fc_k <- tsclust(data_subset, type = "f", k = 2L:5L,
                    preproc = acf_fun, distance = "L2",
                    control = fuzzy_control(version = 1L),
                    seed = 123)
    expect_identical(length(fc_k), 4L)
    fc_k <- lapply(fc_k, reset_nondeterministic)
    assign("fc_k", fc_k, persistent)
})

# ==================================================================================================
# valid input
# ==================================================================================================

test_that("Fuzzy clustering works as expected.", {
    ## ---------------------------------------------------------- univariate fcm
    fcm <- tsclust(data_subset, type = "fuzzy", k = 4L,
                   preproc = acf_fun, distance = "L2",
                   control = fuzzy_control(version = 1L),
                   seed = 123)
    expect_s4_class(fcm, "FuzzyTSClusters")
    fcm <- reset_nondeterministic(fcm)
    assign("fcm", fcm, persistent)

    ## ---------------------------------------------------------- univariate fcmdd
    fcmdd <- tsclust(data_subset, type = "fuzzy", k = 4L,
                     preproc = acf_fun, distance = "L2",
                     centroid = "fcmdd", seed = 123,
                     control = fuzzy_control(version = 1L))
    expect_s4_class(fcmdd, "FuzzyTSClusters")
    fcmdd <- reset_nondeterministic(fcmdd)
    assign("fcmdd", fcmdd, persistent)

    ## ---------------------------------------------------------- multivariate fcm
    dmv <- reinterpolate(data_multivariate, new.length = max(sapply(data_multivariate, NROW)))

    fcm_mv <- tsclust(dmv, type = "fuzzy", k = 4L,
                      distance = "dtw_basic",
                      control = fuzzy_control(version = 1L),
                      seed = 123)
    expect_s4_class(fcm_mv, "FuzzyTSClusters")
    fcm_mv <- reset_nondeterministic(fcm_mv)
    assign("fcm_mv", fcm_mv, persistent)

    ## ---------------------------------------------------------- multivariate fcmdd
    fcmdd_mv <- tsclust(data_multivariate, type = "fuzzy", k = 4L,
                        distance = "dtw_basic", centroid = "fcmdd",
                        control = fuzzy_control(version = 1L),
                        seed = 123)
    expect_s4_class(fcmdd_mv, "FuzzyTSClusters")
    fcmdd_mv <- reset_nondeterministic(fcmdd_mv)
    assign("fcmdd_mv", fcmdd_mv, persistent)
})

# ==================================================================================================
# custom centroid
# ==================================================================================================

test_that("Operations with custom fuzzy centroid complete successfully.", {
    ## ---------------------------------------------------------- with dots
    myfcent <- allcent <- function(x, cl_id, k, cent, cl_old, ...) {
        x <- tslist(x)
        u <- cl_id ^ 2
        cent <- t(u) %*% do.call(rbind, x, TRUE)
        cent <- apply(cent, 2L, "/", e2 = colSums(u))
        tslist(cent)
    }
    expect_output(
        fcent_fcm <- tsclust(data_matrix, k = 20L, type = "fuzzy",
                             distance = "L2", centroid = myfcent,
                             seed = 123, trace = TRUE)
    )
    fcent_fcm <- reset_nondeterministic(fcent_fcm)
    expect_identical(fcent_fcm@centroid, "myfcent")

    ## ---------------------------------------------------------- without dots
    myfcent <- allcent <- function(x, cl_id, k, cent, cl_old) {
        x <- tslist(x)
        u <- cl_id ^ 2
        cent <- t(u) %*% do.call(rbind, x, TRUE)
        cent <- apply(cent, 2L, "/", e2 = colSums(u))
        tslist(cent)
    }
    expect_output(
        fcent_fcm_nd <- tsclust(data_matrix, k = 20L, type = "fuzzy",
                                distance = "L2", centroid = myfcent,
                                seed = 123, trace = TRUE)
    )
    fcent_fcm_nd <- reset_nondeterministic(fcent_fcm_nd)
    expect_identical(fcent_fcm@centroid, "myfcent")

    ## ---------------------------------------------------------- refs
    assign("fcent_fcm", fcent_fcm, persistent)
    assign("fcent_fcm_nd", fcent_fcm_nd, persistent)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
