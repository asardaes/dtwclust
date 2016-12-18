context("\tInvalid inputs")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# invalid data
# =================================================================================================

test_that("Errors in input data are detected", {
    expect_error(dtwclust(NULL), "No data")
    expect_error(dtwclust(NA), "type")
    expect_error(dtwclust(mean), "type")
    expect_error(dtwclust("data"), "type")

    expect_error(dtwclust(as.logical(data_matrix)), "type")

    temp <- data[[1L]]
    temp[2L] <- NA
    data[[1L]] <- temp

    expect_error(dtwclust(data), "missing values")

    data[[1L]] <- numeric()

    expect_error(dtwclust(data), "one point")
})

# =================================================================================================
# invalid distance
# =================================================================================================

test_that("Errors in distance argument are correctly detected.", {
    expect_error(dtwclust(data_matrix, k = 20, distance = mean), "proxy", info = "Function")

    expect_error(dtwclust(data_matrix, k = 20, distance = NULL), "proxy", info = "NULL")

    expect_error(dtwclust(data_matrix, k = 20, distance = NA), "proxy", info = "NA")

    expect_error(dtwclust(data_matrix, k = 20, distance = "dummy"), "proxy", info = "Unregistered")

    expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                 "different length", info = "LBK")

    expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                 "different length", info = "LBI")

    expect_error(dtwclust(data, k = 20, distance = "lbi", control = list(window.size = 18L)),
                 "different length", info = "DTW_LB")
})

# =================================================================================================
# invalid centroid
# =================================================================================================

test_that("Errors in centroid argument are correctly detected.", {
    expect_error(dtwclust(data, centroid = "mean"),
                 "different length", "mean")

    expect_error(dtwclust(data, centroid = "mean"),
                 "different length", "median")

    expect_error(dtwclust(data, centroid = "fcm"),
                 "arg", info = "FCM is for fuzzy")

    expect_error(dtwclust(data, centroid = "goku"),
                 "arg", info = "Unknown centroid")

    expect_error(dtwclust(data, centroid = NULL),
                 "definition", info = "NULL centroid")

    expect_error(dtwclust(data, centroid = NA),
                 "definition", info = "NA centroid")

    expect_error(dtwclust(data, centroid = function(x, ...) { x }),
                 "centroid.*arguments")
})

# =================================================================================================
# invalid control
# =================================================================================================

test_that("Errors in control argument are detected correctly.", {
    expect_error(dtwclust(data, control = NA), "control")
    expect_error(dtwclust(data, control = mean), "control")
    expect_error(dtwclust(data, control = TRUE), "control")
    expect_error(dtwclust(data, control = c(window.size = 15L, iter.max = 15L)), "control")
    expect_error(dtwclust(data, control = list(window.size = c(13L, 14L))),
                 "control", ignore.case = TRUE)
})

# =================================================================================================
# fuzzy clustering
# =================================================================================================

test_that("Invalid combinations in fuzzy clustering are detected.", {
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

    expect_error(dtwclust(data_matrix, type = "f", distance = mean), "proxy", info = "Function")
    expect_error(dtwclust(data_matrix, type = "f", distance = NULL), "proxy", info = "NULL")
    expect_error(dtwclust(data_matrix, type = "f", distance = NA), "proxy", info = "NA")
    expect_error(dtwclust(data_matrix, type = "f", distance = "dummy"), "proxy", info = "Unregistered")

    expect_error(dtwclust(data_subset, type = "f",
                          preproc = reinterpolate, new.length = 205L,
                          distance = "L2", centroid = "mean"))

    expect_error(dtwclust(data_subset, type = "f",
                          preproc = reinterpolate, new.length = 205L,
                          distance = "L2", centroid = "median"))

    expect_error(dtwclust(data_subset, type = "f",
                          preproc = reinterpolate, new.length = 205L,
                          distance = "L2", centroid = "shape"))

    expect_error(dtwclust(data_subset, type = "f",
                          preproc = reinterpolate, new.length = 205L,
                          distance = "L2", centroid = "dba"))

    expect_error(dtwclust(data_subset, type = "f",
                          preproc = reinterpolate, new.length = 205L,
                          distance = "L2", centroid = "pam"))
})

# =================================================================================================
# hierarchical clustering
# =================================================================================================

test_that("Invalid combinations in hierarchical clustering are detected.", {
    expect_error(dtwclust(data, type = "h", k = 101L), "more clusters")

    expect_error(dtwclust(data, type = "h", distance = "lbk"), "different length")
    expect_error(dtwclust(data, type = "h", distance = "lbi"), "different length")
    expect_error(dtwclust(data, type = "h", distance = "dtw_lb"), "different length")

    expect_error(dtwclust(data, type = "h", preproc = "zscore"), "preprocessing")

    expect_error(dtwclust(data_matrix, type = "h", distance = mean), "proxy", info = "Function")
    expect_error(dtwclust(data_matrix, type = "h", distance = NULL), "proxy", info = "NULL")
    expect_error(dtwclust(data_matrix, type = "h", distance = NA), "proxy", info = "NA")
    expect_error(dtwclust(data_matrix, type = "h", distance = "dummy"), "proxy", info = "Unregistered")
})

# =================================================================================================
# partitional clustering
# =================================================================================================

test_that("Invalid combinations in partitional clustering are detected.", {
    expect_error(dtwclust(data, k = 101L), "more clusters")

    expect_error(dtwclust(data, distance = "lbk"), "different length")
    expect_error(dtwclust(data, distance = "lbi"), "different length")
    expect_error(dtwclust(data, distance = "dtw_lb"), "different length")

    expect_error(dtwclust(data, centroid = "mean"), "different length")
    expect_error(dtwclust(data, centroid = "median"), "different length")
    expect_error(dtwclust(data, centroid = "fcm"), "arg")

    expect_error(dtwclust(data, preproc = "zscore"), "preprocessing")
})

# =================================================================================================
# TADPole clustering
# =================================================================================================

test_that("Invalid combinations in tadpole clustering are detected.", {
    expect_error(dtwclust(data, type = "t", k = 20, control = list(window.size = 20L)), "dc")
    expect_error(dtwclust(data, type = "t", k = 20, dc = 1.5), "window.size")
    expect_error(dtwclust(data, type = "t", k = 20, dc = 1.5, control = list(window.size = 20L)),
                 "same length")
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
