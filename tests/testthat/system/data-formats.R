context("\tData formats")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# data formats
# =================================================================================================

test_that("Results are the same regardless of data format, as long as supported.", {
    pc_matrix <- dtwclust(data_matrix[1L:20L, ], type = "partitional", k = 4L,
                          distance = "L2", centroid = "pam",
                          seed = 939)

    pc_matrix <- reset_nondeterministic(pc_matrix)

    pc_list <- dtwclust(data_reinterpolated_subset, type = "partitional", k = 4L,
                        distance = "L2", centroid = "pam",
                        seed = 939)

    pc_list <- reset_nondeterministic(pc_list)

    pc_df<- dtwclust(as.data.frame(data_matrix[1L:20L, ]),
                     type = "partitional", k = 4L,
                     distance = "L2", centroid = "pam",
                     seed = 939)

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

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
