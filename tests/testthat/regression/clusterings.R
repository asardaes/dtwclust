context("\tClusterings")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# fuzzy
# =================================================================================================

with(persistent, {
    test_that("Fuzzy clustering gives the same results as references.", {
        skip_on_cran()

        expect_equal_to_reference(fc_k, file_name(fc_k))
        expect_equal_to_reference(fcm, file_name(fcm))
        expect_equal_to_reference(fcmdd, file_name(fcmdd))
        expect_equal_to_reference(fcm_mv, file_name(fcm_mv))
        expect_equal_to_reference(fcmdd_mv, file_name(fcmdd_mv))
    })
})

# =================================================================================================
# hierarchical
# =================================================================================================

with(persistent, {
    test_that("Hierarchical clustering gives the same results as references.", {
        skip_on_cran()

        expect_equal_to_reference(hc_k, file_name(hc_k))
        expect_equal_to_reference(hc_all, file_name(hc_all, x32 = TRUE))
        expect_equal_to_reference(hc_lbi, file_name(hc_lbi))
        expect_equal_to_reference(hc_cent, file_name(hc_cent))
        expect_equal_to_reference(hc_diana, file_name(hc_diana, x32 = TRUE))
    })
})

# =================================================================================================
# partitional
# =================================================================================================

with(persistent, {
    test_that("Partitional clustering gives the same results as references.", {
        skip_on_cran()

        expect_equal_to_reference(pc_k, file_name(pc_k))
        expect_equal_to_reference(pc_rep, file_name(pc_rep))
        expect_equal_to_reference(pc_krep, file_name(pc_krep))

        expect_equal_to_reference(pc_dtwb, file_name(pc_dtwb))
        expect_equal_to_reference(pc_dtwb_npampre, file_name(pc_dtwb_npampre))
        expect_equal_to_reference(pc_dtwb_distmat, file_name(pc_dtwb_distmat))
        expect_equal_to_reference(pc_dtwlb, file_name(pc_dtwlb))

        expect_equal_to_reference(pc_kshape, file_name(pc_kshape))
        expect_equal_to_reference(pc_dba, file_name(pc_dba))
        expect_equal_to_reference(pc_mv_pam, file_name(pc_mv_pam))
        expect_equal_to_reference(pc_mv_dba, file_name(pc_mv_dba, x32 = TRUE))

        expect_equal_to_reference(pc_tadp, file_name(pc_tadp))
        expect_equal_to_reference(pc_tadp_lbi, file_name(pc_tadp_lbi))
        expect_equal_to_reference(pc_tadp_cent, file_name(pc_tadp_cent))

        expect_equal_to_reference(pc_cr, file_name(pc_cr))
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
