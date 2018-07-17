context("    Invalid inputs")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# invalid data
# ==================================================================================================

test_that("Errors in input data are detected by tsclust.", {
    expect_error(tsclust(NULL), "data")
    expect_error(tsclust(NA), "type")
    expect_error(tsclust(mean), "type")
    expect_error(tsclust("data"), "type")

    expect_error(tsclust(as.logical(data_matrix)), "type")

    temp <- data[[1L]]
    temp[2L] <- NA
    data[[1L]] <- temp

    expect_error(tsclust(data), "missing values")

    data[[1L]] <- numeric()

    expect_error(tsclust(data), "one point")
})

# ==================================================================================================
# invalid control
# ==================================================================================================

test_that("Errors in control argument are detected correctly by tsclust.", {
    expect_error(tsclust(data, control = NA), "control")
    expect_error(tsclust(data, control = mean), "control")
    expect_error(tsclust(data, control = TRUE), "control")
    expect_error(tsclust(data, control = partitional_control(distmat = matrix(0, 2L, 2L))))
})

# ==================================================================================================
# fuzzy clustering
# ==================================================================================================

test_that("Invalid combinations in fuzzy clustering are detected by tsclust.", {
    expect_error(tsclust(data_matrix, type = "f", k = 101L), "more clusters")

    args <- tsclust_args(dist = list(window.size = 15L))

    expect_error(tsclust(data, type = "f", distance = "lbk", args = args), "different length")
    expect_error(tsclust(data, type = "f", distance = "lbi", args = args), "different length")
    expect_error(tsclust(data, type = "f", distance = "dtw_lb", args = args), "different length")
    expect_error(tsclust(data, type = "f", distance = "dtw"), "different length")
    expect_error(tsclust(data, type = "f", distance = "dtw2"), "different length")
    expect_error(tsclust(data, type = "f", distance = "dtw_basic"), "different length")
    expect_error(tsclust(data, type = "f", distance = "sbd"), "different length")

    expect_error(tsclust(data_matrix, type = "f", preproc = "zscore"), "preprocessing")

    expect_error(tsclust(data_matrix, type = "f", distance = mean), "proxy", info = "Function")
    expect_error(tsclust(data_matrix, type = "f", distance = NULL), "proxy", info = "NULL")
    expect_error(tsclust(data_matrix, type = "f", distance = NA), "proxy", info = "NA")
    expect_error(tsclust(data_matrix, type = "f", distance = "dummy"), "proxy", info = "Unregistered")

    expect_error(tsclust(data_subset, type = "f",
                         preproc = reinterpolate, new.length = 205L,
                         distance = "L2", centroid = "mean"))

    expect_error(tsclust(data_subset, type = "f",
                         preproc = reinterpolate, new.length = 205L,
                         distance = "L2", centroid = "median"))

    expect_error(tsclust(data_subset, type = "f",
                         preproc = reinterpolate, new.length = 205L,
                         distance = "L2", centroid = "shape"))

    expect_error(tsclust(data_subset, type = "f",
                         preproc = reinterpolate, new.length = 205L,
                         distance = "L2", centroid = "dba"))

    expect_error(tsclust(data_subset, type = "f",
                         preproc = reinterpolate, new.length = 205L,
                         distance = "L2", centroid = "pam"))
})

# ==================================================================================================
# hierarchical clustering
# ==================================================================================================

test_that("Invalid combinations in hierarchical clustering are detected by tsclust.", {
    expect_error(tsclust(data, type = "h", k = 101L), "more clusters")

    expect_error(tsclust(data, type = "h", distance = "lbk"), "different length")
    expect_error(tsclust(data, type = "h", distance = "lbi"), "different length")
    expect_error(tsclust(data, type = "h", distance = "dtw_lb"), "different length")

    expect_error(tsclust(data, type = "h", preproc = "zscore"), "preprocessing")

    expect_error(tsclust(data_matrix, type = "h", distance = mean), "proxy", info = "Function")
    expect_error(tsclust(data_matrix, type = "h", distance = NULL), "proxy", info = "NULL")
    expect_error(tsclust(data_matrix, type = "h", distance = NA), "proxy", info = "NA")
    expect_error(tsclust(data_matrix, type = "h", distance = "dummy"), "proxy", info = "Unregistered")
})

# ==================================================================================================
# partitional clustering
# ==================================================================================================

test_that("Invalid combinations in partitional clustering are detected by tsclust.", {
    expect_error(tsclust(data, k = 101L), "more clusters")

    expect_error(tsclust(data_matrix, k = 20, distance = mean), "proxy", info = "Function")
    expect_error(tsclust(data_matrix, k = 20, distance = NULL), "proxy", info = "NULL")
    expect_error(tsclust(data_matrix, k = 20, distance = NA), "proxy", info = "NA")
    expect_error(tsclust(data_matrix, k = 20, distance = "dummy"), "proxy", info = "Unregistered")
    expect_error(tsclust(data, k = 20, distance = "lbk",
                         args = tsclust_args(dist = list(window.size = 18L))),
                 "different length", info = "LBK")
    expect_error(tsclust(data, k = 20, distance = "lbi",
                         args = tsclust_args(dist = list(window.size = 18L))),
                 "different length", info = "LBK")
    expect_error(tsclust(data, k = 20, distance = "dtw_lb",
                         args = tsclust_args(dist = list(window.size = 18L))),
                 "different length", info = "DTW_LB")

    expect_error(tsclust(data, centroid = "mean"),
                 "different length", info = "mean")
    expect_error(tsclust(data, centroid = "median"),
                 "different length", info = "median")
    expect_error(tsclust(data, centroid = "fcm"),
                 "arg", info = "FCM is for fuzzy")
    expect_error(tsclust(data, centroid = "fcmdd"),
                 "arg", info = "FCMdd is for fuzzy")
    expect_error(tsclust(data, centroid = "goku"),
                 "arg", info = "Unknown centroid")
    expect_error(tsclust(data, centroid = NULL),
                 "definition", info = "NULL centroid")
    expect_error(tsclust(data, centroid = NA),
                 "definition", info = "NA centroid")
    expect_error(tsclust(data, centroid = function(x, ...) { x }),
                 "centroid.*arguments",
                 info = "custom centroid function requires specific formal parameters")

    expect_error(tsclust(data, preproc = "zscore"), "preprocessing",
                 info = "got char instead of function")
})

# ==================================================================================================
# TADPole clustering
# ==================================================================================================

test_that("Invalid combinations in tadpole clustering are detected by tsclust.", {
    expect_error(tsclust(data, type = "t", k = 20, control = tadpole_control(window.size = 20L)),
                 "dc")
    expect_error(tsclust(data, type = "t", k = 20, control = tadpole_control(dc = 1.5)),
                 "window.size")
    expect_error(tsclust(data, type = "t", k = 20,
                         control = tadpole_control(dc = 1.5, window.size = 20L)),
                 "same length")
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
