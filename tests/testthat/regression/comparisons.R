context("\tCompare clusterings")

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

        expect_equal_to_reference(comp_all, file_name(comp_all))
        expect_equal_to_reference(comp_gak, file_name(comp_gak))
        expect_equal_to_reference(comp_dba, file_name(comp_dba))
    })
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
