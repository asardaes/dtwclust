context("\tCVIs")

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

        expect_equal_to_reference(base_cvis, file_name(base_cvis))
        expect_equal_to_reference(internal_fcvis, file_name(internal_fcvis))
        expect_equal_to_reference(external_fcvis, file_name(external_fcvis))
        expect_equal_to_reference(cvis_tadp, file_name(cvis_tadp))
        expect_equal_to_reference(cvis_hc, file_name(cvis_hc))
        expect_equal_to_reference(cvis_tadp_cent, file_name(cvis_tadp_cent))
        expect_equal_to_reference(cvis_hc_cent, file_name(cvis_hc_cent))
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
