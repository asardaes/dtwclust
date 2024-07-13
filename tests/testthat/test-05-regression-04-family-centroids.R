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
        local_edition(2)

        expect_known_value(cent_mean, file_name(cent_mean), info = "Univariate")
        expect_known_value(cent_mv_mean, file_name(cent_mv_mean), info = "Multivariate")
        expect_known_value(cent_median, file_name(cent_median), info = "Univariate")
        expect_known_value(cent_mv_median, file_name(cent_mv_median), info = "Multivariate")
        expect_known_value(cent_shape, file_name(cent_shape), info = "Univariate")
        expect_known_value(cent_mv_shape, file_name(cent_mv_shape), info = "Multivariate")
        expect_known_value(cent_pam, file_name(cent_pam), info = "Univariate without distmat")
        expect_known_value(cent_mv_pam, file_name(cent_mv_pam), info = "Multivariate")

        tol <- if (.Machine$sizeof.pointer == 4L) 1e-2 else sqrt(.Machine$double.eps)
        expect_known_value(cent_dba, file_name(cent_dba, x32 = TRUE), tolerance = tol, info = "Univariate")
        expect_known_value(cent_mv_dba, file_name(cent_mv_dba, x32 = TRUE), tolerance = tol, info = "Multivariate")

        expect_known_value(cent_mv_dba_bys, file_name(cent_mv_dba_bys), info = "DBA by series")
    })

    test_that("Centroids calculated with SDTWC families give the same results as references.", {
        skip_on_cran()
        skip_if(tolower(Sys.info()[["sysname"]]) == "windows" & isTRUE(as.logical(Sys.getenv("CI"))), "On Windows CI")
        local_edition(2)

        expect_known_value(cent_sdtwc, file_name(cent_sdtwc), tolerance = 1e-6, info = "SDTWC Univariate")
        expect_known_value(cent_mv_sdtwc, file_name(cent_mv_sdtwc), info = "SDTWC Multivariate")
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
