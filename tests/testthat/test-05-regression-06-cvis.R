# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# centroids
# =================================================================================================

with(persistent, {
    test_that("CVIs give the same results as references.", {
        skip_on_cran()
        local_edition(2)

        expect_known_value(base_cvis, file_name(base_cvis))
        expect_known_value(internal_fcvis, file_name(internal_fcvis))
        expect_known_value(external_fcvis, file_name(external_fcvis))
        expect_known_value(cvis_tadp, file_name(cvis_tadp))
        expect_known_value(cvis_hc, file_name(cvis_hc))
        expect_known_value(cvis_tadp_cent, file_name(cvis_tadp_cent))
        expect_known_value(cvis_hc_cent, file_name(cvis_hc_cent))
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
