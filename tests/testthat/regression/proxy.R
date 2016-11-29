context("\tProxy distances")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

## Distances to test
included_distances <- c("lbk", "lbi", "sbd", "dtw_basic", "dtw_lb", "gak")

## Data
x <- data_reinterpolated[3L:8L]

# =================================================================================================
# pairwise
# =================================================================================================

test_that("Pairwise proxy distances give the same result as references", {
    for (distance in included_distances) {
        d <- proxy::dist(x[1L:3L], x[4L:6L], method = distance,
                         window.size = 15L, sigma = 100,
                         pairwise = TRUE)

        expect_equal_to_reference(d, paste0("rds/pdist_", distance, ".rds"))
    }
})

# =================================================================================================
# multivariate
# =================================================================================================

test_that("Included (valid) distances can accept multivariate series.", {
    skip_on_cran()

    for (distance in c("dtw_basic", "gak")) {
        mv <- proxy::dist(data_multivariate, method = distance,
                          window.size = 18L, sigma = 100)

        expect_equal_to_reference(mv, paste0("rds/mv_", distance, ".rds"))
    }
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
