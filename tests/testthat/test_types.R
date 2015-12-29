context("Test data types")

# =================================================================================================
# matrix
# =================================================================================================

pc_matrix <- dtwclust(data_matrix, type = "partitional", k = 20,
                      distance = "L2", centroid = "pam",
                      preproc = NULL, seed = 123,
                      control = ctrl)

pc_matrix <- reset_nondeterministic(pc_matrix)

test_that("matrix input gives the same result as reference",
          my_expect_equal_to_reference(pc_matrix))

# =================================================================================================
# list
# =================================================================================================

pc_list <- dtwclust(data_list, type = "partitional", k = 20,
                    distance = "L1", centroid = "pam",
                    preproc = NULL, seed = 123,
                    conrol = ctrl)

pc_list <- reset_nondeterministic(pc_list)

test_that("list input gives the same result as reference",
          my_expect_equal_to_reference(pc_list))

# =================================================================================================
# data.frame
# =================================================================================================

test_that("logical input gives data type error",
          expect_error(dtwclust(as.logical(data_matrix), k = 20, distance = "L2"),
                       "type"))
