context("    Clusterings")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# fuzzy
# ==================================================================================================

with(persistent, {
    test_that("Fuzzy clustering gives the same results as references.", {
        skip_on_cran()

        expect_known_value(fc_k, file_name(fc_k))
        expect_known_value(fcm, file_name(fcm))
        expect_known_value(fcmdd, file_name(fcmdd))
        expect_known_value(fcm_mv, file_name(fcm_mv))
        expect_known_value(fcmdd_mv, file_name(fcmdd_mv))

        # notice files are the same, results should be equal
        expect_known_value(fcent_fcm, file_name(fcent_fcm), info = "Custom fuzzy c-means")
        expect_known_value(fcent_fcm_nd, file_name(fcent_fcm), info = "Custom fuzzy c-means")
    })
})

# ==================================================================================================
# hierarchical
# ==================================================================================================

with(persistent, {
    test_that("Hierarchical clustering gives the same results as references.", {
        skip_on_cran()

        expect_known_value(hc_k, file_name(hc_k))
        expect_known_value(hc_all, file_name(hc_all))
        expect_known_value(hc_lbi, file_name(hc_lbi))
        expect_known_value(hc_cent, file_name(hc_cent))
        expect_known_value(hc_cent2, file_name(hc_cent2))
        expect_known_value(hc_diana, file_name(hc_diana))
    })
})

# ==================================================================================================
# partitional
# ==================================================================================================

with(persistent, {
    test_that("Partitional clustering gives the same results as references.", {
        skip_on_cran()

        expect_known_value(pc_k, file_name(pc_k))
        expect_known_value(pc_rep, file_name(pc_rep))
        expect_known_value(pc_krep, file_name(pc_krep))

        expect_known_value(pc_dtwb, file_name(pc_dtwb))
        expect_known_value(pc_dtwb_npampre, file_name(pc_dtwb_npampre))
        expect_known_value(pc_dtwb_distmat, file_name(pc_dtwb_distmat))
        expect_known_value(pc_dtwlb, file_name(pc_dtwlb))

        expect_known_value(pc_kshape, file_name(pc_kshape))
        expect_known_value(pc_dba, file_name(pc_dba))
        expect_known_value(pc_mv_pam, file_name(pc_mv_pam))
        expect_known_value(pc_mv_dba, file_name(pc_mv_dba))
        expect_known_value(pc_sdtw, file_name(pc_sdtw))

        expect_known_value(pc_tadp, file_name(pc_tadp))
        expect_known_value(pc_tadp_lbi, file_name(pc_tadp_lbi))
        expect_known_value(pc_tadp_cent, file_name(pc_tadp_cent))

        expect_known_value(pc_cr, file_name(pc_cr))

        # notice files are the same, results should be equal
        expect_known_value(cent_colMeans, file_name(cent_colMeans), info = "Custom colMeans")
        expect_known_value(cent_colMeans_nd, file_name(cent_colMeans), info = "Custom colMeans")
    })
})

# ==================================================================================================
# clean
# ==================================================================================================
rm(list = setdiff(ls(), ols))
