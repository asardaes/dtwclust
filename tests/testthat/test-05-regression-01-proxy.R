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
    skip_on_cran()
    local_edition(2)
    for (distance in included_distances) {
        d <- proxy::dist(x[1L:3L], x[4L:6L], method = distance,
                         window.size = 15L, sigma = 100,
                         pairwise = TRUE)

        expect_known_value(d, paste0("rds/pdist_", distance, ".rds"),
                           check.attributes = FALSE)
    }
})

# =================================================================================================
# multivariate
# =================================================================================================

test_that("Included (valid) distances can accept multivariate series.", {
    skip_on_cran()
    local_edition(2)

    for (distance in c("dtw_basic", "gak")) {
        mv <- proxy::dist(data_multivariate, method = distance,
                          window.size = 18L, sigma = 100)

        expect_known_value(mv, paste0("rds/mv_", distance, ".rds"))
    }
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
