# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# distmats
# =================================================================================================

with(persistent, {
    test_that("Distance matrices calculated with families give the same results as references.", {
        skip_on_cran()
        local_edition(2)

        expect_known_value(distmat_lbk, file_name(distmat_lbk), info = "LBK")
        expect_known_value(distmat_lbi, file_name(distmat_lbi), info = "LBI")
        expect_known_value(distmat_sbd, file_name(distmat_sbd), info = "SBD")
        expect_known_value(distmat_dtwlb, file_name(distmat_dtwlb), info = "DTW_LB")
        expect_known_value(distmat_dtw, file_name(distmat_dtw), info = "DTW")
        expect_known_value(distmat_dtw2, file_name(distmat_dtw2), info = "DTW2")
        expect_known_value(distmat_dtwb, file_name(distmat_dtwb), info = "DTW_BASIC")
        expect_known_value(distmat_gak, file_name(distmat_gak), info = "GAK")
        expect_known_value(distmat_sdtw, file_name(distmat_sdtw), info = "SDTW")
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
