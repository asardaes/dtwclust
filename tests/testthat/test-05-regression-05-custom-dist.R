# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# centroids
# =================================================================================================

with(persistent, {
    test_that("A custom distance in dtwclust give the same results as references.", {
        skip_on_cran()
        local_edition(2)

        expect_known_value(pc_ndtw, file_name(pc_ndtw), info = "nDTW")
        expect_known_value(pc_ndtw_sym, file_name(pc_ndtw_sym), info = "Symmetric nDTW")
        expect_known_value(pc_ndtw_par, file_name(pc_ndtw_par), info = "Params with nDTW")
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
