context("\tDistances")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

## Functions
functions <- c("lb_keogh", "lb_improved", "dtw_lb", "SBD", "dtw_basic", "GAK")

## Support
supports_mv <- c("dtw_basic", "GAK")
supports_diff_lengths <- c("SBD", "dtw_basic", "GAK")

## Extra arguments
args <- list(
    lb_keogh = list(
        list(window.size = 15L, norm = "L1", force.symmetry = FALSE),
        list(window.size = 15L, norm = "L1", force.symmetry = TRUE),
        list(window.size = 15L, norm = "L2", force.symmetry = FALSE),
        list(window.size = 15L, norm = "L2", force.symmetry = TRUE)
    ),
    lb_improved = list(
        list(window.size = 15L, norm = "L1", force.symmetry = FALSE),
        list(window.size = 15L, norm = "L1", force.symmetry = TRUE),
        list(window.size = 15L, norm = "L2", force.symmetry = FALSE),
        list(window.size = 15L, norm = "L2", force.symmetry = TRUE)
    ),
    dtw_lb = list(
        list(window.size = 15L, norm = "L1"),
        list(window.size = 15L, norm = "L2")
    ),
    SBD = list(
    ),
    dtw_basic = list(
        list(window.size = 15L, norm = "L1", step.pattern = symmetric1),
        list(window.size = 15L, norm = "L1", step.pattern = symmetric2),
        list(window.size = 15L, norm = "L2", step.pattern = symmetric1),
        list(window.size = 15L, norm = "L2", step.pattern = symmetric2),
        list(window.size = 15L, norm = "L1", step.pattern = symmetric2, normalize = TRUE),
        list(window.size = 15L, norm = "L2", step.pattern = symmetric2, normalize = TRUE)
    ),
    GAK = list(
        sigma = 100, window.size = 15L,
        sigma = 100, window.size = NULL,
        sigma = 100, normalize = FALSE)
)

## Univariate series
x_uv <- data[[1L]]
y_uv_same_length <- data[[2L]]
y_uv_diff_length <- data[[100L]]

## Multivariate series
x_mv <- data_multivariate[[1L]]
y_mv_same_length <- data_multivariate[[2L]]
y_mv_diff_length <- data_multivariate[[20L]]

## Invalid
invalid_inputs <- list(
    "NA" = NA,
    "NULL" = NULL,
    "char" = "1",
    "bool" = TRUE,
    "empty" = numeric(),
    "miss" = c(1, NA)
)

# =================================================================================================
# invalid inputs
# =================================================================================================

test_that("Invalid inputs are detected correctly in the distance functions.", {
    for (foo in functions) {
        for (input in names(invalid_inputs)) {
            expect_error(do.call(foo, list(x = invalid_inputs[[input]], y = x_uv)),
                         info = paste(foo, "=>", input, "in x"))

            expect_error(do.call(foo, list(x = x_uv, y = invalid_inputs[[input]])),
                         info = paste(foo, "=>", input, "in y"))

            if (!(foo %in% supports_mv)) {
                expect_error(do.call(foo, list(x = x_mv, y = x_mv)),
                             info = paste(foo, "=> multivariate"))
            }

            if (!(foo %in% supports_diff_lengths)) {
                expect_error(do.call(foo, list(x = x_uv, y = y_uv_diff_length)),
                             info = paste(foo, "=> different lengths"))
            }
        }
    }
})

# =================================================================================================
# valid inputs
# =================================================================================================

test_that("Valid inputs provide a result greater than zero", {
    for (foo in functions) {
        for (arg in args[[foo]]) {
            distance_value <- do.call(foo,
                                      c(list(x = x_uv,
                                             y = y_uv_same_length),
                                        arg))

            if (foo == "lb_keogh") distance_value <- distance_value$d

            expect_gt(distance_value, 0, label = paste0("distance with ", foo))

            if (mv <- foo %in% supports_mv) {
                distance_value <- do.call(foo,
                                          c(list(x = x_mv,
                                                 y = y_mv_same_length),
                                            arg))

                expect_gt(distance_value, 0, label = paste0("multivariate distance with ", foo))
            }

            if (dl <- foo %in% supports_diff_lengths) {
                distance_value <- do.call(foo,
                                          c(list(x = x_uv,
                                                 y = y_uv_diff_length),
                                            arg))

                expect_gt(distance_value, 0, label = paste0("distance with different lengths with ", foo))
            }

            if (mv && dl) {
                distance_value <- do.call(foo,
                                          c(list(x = x_mv,
                                                 y = y_mv_diff_length),
                                            arg))

                expect_gt(distance_value, 0, label = paste0("multivariate distance with different lengths with ", foo))
            }
        }
    }
})

# =================================================================================================
# GAK's sigma estimation
# =================================================================================================

test_that("GAK can estimate sigma.", {
    gak_distance <- GAK(data[[1L]], data[[100L]])
    expect_gt(attr(gak_distance, "sigma"), 0)
})

# =================================================================================================
# dtw_lb with dtw::dtw
# =================================================================================================

test_that("dtw_lb has the same result regardless of dtw.func.", {
    distmat_with_dtwbasic <- dtw_lb(data_reinterpolated[1L:50L], data_reinterpolated[51L:100L],
                                    window.size = 15L, step.pattern = symmetric1)
    distmat_with_dtw <- dtw_lb(data_reinterpolated[1L:50L], data_reinterpolated[51L:100L],
                               window.size = 15L, step.pattern = symmetric1, dtw.func = "dtw")
    expect_equal(distmat_with_dtwbasic, distmat_with_dtw)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
