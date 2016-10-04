context("Test data types")

# =================================================================================================
# data errors
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
# data formats
# =================================================================================================

test_that("Results are the same regardless of data format, as long as supported.", {
    pc_matrix <- dtwclust(data_matrix[1L:20L, ], type = "partitional", k = 4L,
                          distance = "L2", centroid = "pam",
                          seed = 939)

    pc_matrix <- reset_nondeterministic(pc_matrix)

    pc_list <- dtwclust(data_reinterpolated[1L:20L], type = "partitional", k = 4L,
                        distance = "L2", centroid = "pam",
                        seed = 939)

    pc_list <- reset_nondeterministic(pc_list)

    suppressMessages(
        pc_df<- dtwclust(as.data.frame(data_matrix[1L:20L, ]),
                         type = "partitional", k = 4L,
                         distance = "L2", centroid = "pam",
                         seed = 939)
    )

    pc_df <- reset_nondeterministic(pc_df)

    pc_matrix <- lapply(pc_matrix, function(obj) {
        obj@call <- as.call(list("zas", a = 1))
    })

    pc_list <- lapply(pc_list, function(obj) {
        obj@call <- as.call(list("zas", a = 1))
    })

    pc_df <- lapply(pc_df, function(obj) {
        obj@call <- as.call(list("zas", a = 1))
    })

    expect_identical(pc_matrix, pc_list)
    expect_identical(pc_matrix, pc_df)
})
