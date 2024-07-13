# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# comparisons
# =================================================================================================

with(persistent, {
    test_that("Compare clusterings gives the same results as references.", {
        skip_on_cran()
        local_edition(2)

        expect_known_value(comp_all, file_name(comp_all))
        expect_known_value(comp_gak, file_name(comp_gak))
        expect_known_value(comp_dba, file_name(comp_dba))
        expect_known_value(comp_sdtwc, file_name(comp_sdtwc))
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
