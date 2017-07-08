context("\tFamilies' centroids")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# centroids
# =================================================================================================

with(persistent, {
    test_that("Centroids calculated with families give the same results as references.", {
        skip_on_cran()

        expect_equal_to_reference(cent_mean, file_name(cent_mean), info = "Univariate")
        expect_equal_to_reference(cent_mv_mean, file_name(cent_mv_mean), info = "Multivariate")
        expect_equal_to_reference(cent_median, file_name(cent_median), info = "Univariate")
        expect_equal_to_reference(cent_mv_median, file_name(cent_mv_median), info = "Multivariate")
        expect_equal_to_reference(cent_shape, file_name(cent_shape), info = "Univariate")
        expect_equal_to_reference(cent_mv_shape, file_name(cent_mv_shape), info = "Multivariate")
        expect_equal_to_reference(cent_pam, file_name(cent_pam), info = "Univariate without distmat")
        expect_equal_to_reference(cent_mv_pam, file_name(cent_mv_pam), info = "Multivariate")
        expect_equal_to_reference(cent_dba, file_name(cent_dba, x32 = TRUE), info = "Univariate")
        expect_equal_to_reference(cent_mv_dba, file_name(cent_mv_dba, x32 = TRUE), info = "Multivariate")
        expect_equal_to_reference(cent_mv_dba_bys, file_name(cent_mv_dba_bys), info = "DBA by series")

        ## notice files are the same, results should be equal
        expect_equal_to_reference(cent_colMeans, file_name(cent_colMeans), info = "Custom colMeans")
        expect_equal_to_reference(cent_colMeans_nd, file_name(cent_colMeans), info = "Custom colMeans")
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
