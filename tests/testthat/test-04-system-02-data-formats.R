# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# data formats
# ==================================================================================================

test_that("Results are the same regardless of data format, as long as supported by tsclust.", {
    pc_matrix <- tsclust(data_matrix[1L:20L, ], type = "partitional", k = 4L,
                         distance = "L2", centroid = "pam",
                         seed = 939)

    pc_matrix <- reset_nondeterministic(pc_matrix)

    pc_list <- tsclust(data_reinterpolated_subset, type = "partitional", k = 4L,
                       distance = "L2", centroid = "pam",
                       seed = 939)

    pc_list <- reset_nondeterministic(pc_list)

    df <- as.data.frame(data_matrix[1L:20L, ])
    colnames(df) <- NULL
    pc_df<- tsclust(df, type = "partitional", k = 4L,
                    distance = "L2", centroid = "pam",
                    seed = 939)

    pc_df <- reset_nondeterministic(pc_df)

    pc_matrix@call <- pc_list@call <- pc_df@call <- as.call(list("zas", a = 1))

    expect_identical(pc_matrix, pc_list)
    expect_identical(pc_matrix, pc_df)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
