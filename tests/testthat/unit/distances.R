context("    Included distances")

# ==================================================================================================
# setup
# ==================================================================================================

# Original objects in env
ols <- ls()

# Functions
functions <- c("lb_keogh", "lb_improved", "dtw_lb", "SBD", "dtw_basic", "GAK", "sdtw")

# Support
supports_mv <- c("dtw_basic", "GAK", "sdtw")
supports_diff_lengths <- c("SBD", "dtw_basic", "GAK", "sdtw")

# Extra arguments
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
        list(znorm = FALSE),
        list(znorm = TRUE)
    ),
    dtw_basic = list(
        list(window.size = 15L, norm = "L1", step.pattern = dtw::symmetric1),
        list(window.size = 15L, norm = "L1", step.pattern = dtw::symmetric2),
        list(window.size = 15L, norm = "L2", step.pattern = dtw::symmetric1),
        list(window.size = 15L, norm = "L2", step.pattern = dtw::symmetric2),
        list(window.size = 15L, norm = "L1", step.pattern = dtw::symmetric2, normalize = TRUE),
        list(window.size = 15L, norm = "L2", step.pattern = dtw::symmetric2, normalize = TRUE)
    ),
    GAK = list(
        list(sigma = 100, window.size = 15L),
        list(sigma = 100, window.size = NULL),
        list(sigma = 100, normalize = FALSE)
    ),
    sdtw = list(
        list()
    )
)

# Univariate series
x_uv <- data[[1L]]
y_uv_same_length <- data[[2L]]
y_uv_diff_length <- data[[100L]]

# Multivariate series
x_mv <- data_multivariate[[1L]]
y_mv_same_length <- data_multivariate[[2L]]
y_mv_diff_length <- data_multivariate[[20L]]

# Invalid
invalid_inputs <- list(
    "NA" = NA,
    "NULL" = NULL,
    "char" = "1",
    "bool" = TRUE,
    "empty" = numeric(),
    "miss" = c(1, NA)
)

# ==================================================================================================
# invalid inputs
# ==================================================================================================

test_that("Invalid inputs are detected correctly in the distance functions.", {
    for (foo in functions) {
        for (input in names(invalid_inputs)) {
            expect_error(do.call(foo, list(x = invalid_inputs[[input]], y = x_uv)),
                         info = paste("function", foo, "with", input, "input in x"))

            expect_error(do.call(foo, list(x = x_uv, y = invalid_inputs[[input]])),
                         info = paste("function", foo, "with", input, "input in y"))

            if (!(foo %in% supports_mv)) {
                expect_error(do.call(foo, list(x = x_mv, y = x_mv)),
                             info = paste("function", foo, "with multivariate input"))
            }

            if (!(foo %in% supports_diff_lengths)) {
                expect_error(do.call(foo, list(x = x_uv, y = y_uv_diff_length)),
                             info = paste("function", foo, "with different-length input"))
            }
        }
    }
})

# ==================================================================================================
# valid inputs
# ==================================================================================================

test_that("Valid inputs provide a result different than zero", {
    for (foo in functions) {
        for (arg in args[[foo]]) {
            distance_value <- do.call(foo,
                                      c(list(x = x_uv,
                                             y = y_uv_same_length),
                                        arg))

            if (foo %in% c("lb_keogh", "SBD")) distance_value <- distance_value$d
            expect_true(distance_value != 0, label = paste0("distance with ", foo))

            if (mv <- foo %in% supports_mv) {
                distance_value <- do.call(foo,
                                          c(list(x = x_mv,
                                                 y = y_mv_same_length),
                                            arg))

                expect_true(distance_value != 0, label = paste0("multivariate distance with ", foo))
            }

            if (dl <- foo %in% supports_diff_lengths) {
                distance_value <- do.call(foo,
                                          c(list(x = x_uv,
                                                 y = y_uv_diff_length),
                                            arg))

                if (foo %in% c("SBD")) distance_value <- distance_value$d
                expect_true(distance_value != 0, label = paste0("distance with different lengths with ", foo))
            }

            if (mv && dl) {
                distance_value <- do.call(foo,
                                          c(list(x = x_mv,
                                                 y = y_mv_diff_length),
                                            arg))

                expect_true(distance_value != 0, label = paste0("multivariate distance with different lengths with ", foo))
            }
        }
    }
})

# ==================================================================================================
# DTW's lower bounds
# ==================================================================================================

test_that("Passing pre-computed envelopes to the lower bounds works correctly.", {
    envelopes <- compute_envelope(y_uv_diff_length, window.size = 15L)
    wrong_lower_envelope <- envelopes$lower
    wrong_upper_envelope <- envelopes$upper

    envelopes <- compute_envelope(y_uv_same_length, window.size = 15L)
    correct_lower_envelope <- envelopes$lower
    correct_upper_envelope <- envelopes$upper

    expect_error(lb_keogh(x_uv,
                          lower.env = wrong_lower_envelope,
                          upper.env = correct_upper_envelope),
                 regexp = "mismatch.*envelope")
    expect_error(lb_keogh(x_uv,
                          lower.env = correct_lower_envelope,
                          upper.env = wrong_upper_envelope),
                 regexp = "mismatch.*envelope")
    expect_gt(lb_keogh(x_uv,
                       lower.env = correct_lower_envelope,
                       upper.env = correct_upper_envelope)$d,
              0)

    expect_error(lb_improved(x_uv, y_uv_same_length, window.size = 15L,
                             lower.env = wrong_lower_envelope,
                             upper.env = correct_upper_envelope),
                 regexp = "mismatch.*envelope")
    expect_error(lb_improved(x_uv, y_uv_same_length, window.size = 15L,
                             lower.env = correct_lower_envelope,
                             upper.env = wrong_upper_envelope),
                 regexp = "mismatch.*envelope")
    expect_gt(lb_improved(x_uv, y_uv_same_length, window.size = 15L,
                          lower.env = correct_lower_envelope,
                          upper.env = correct_upper_envelope),
              0)
})

# ==================================================================================================
# dtw_lb
# ==================================================================================================

test_that("dtw_lb gives the same result regardless of dtw.func.", {
    distmat_with_dtwbasic <- dtw_lb(data_reinterpolated[1L:50L], data_reinterpolated[51L:100L],
                                    window.size = 15L, step.pattern = dtw::symmetric1)
    distmat_with_dtw <- dtw_lb(data_reinterpolated[1L:50L], data_reinterpolated[51L:100L],
                               window.size = 15L, step.pattern = dtw::symmetric1, dtw.func = "dtw")
    expect_equal(distmat_with_dtwbasic, distmat_with_dtw, check.attributes = FALSE)
})

test_that("dtw_lb gives the same result for different nn.margin and corresponding inputs.", {
    # Calculate the DTW distance between a certain subset aided with the lower bound
    d <- dtw_lb(data_reinterpolated[1:5], data_reinterpolated[6:50],
                window.size = 20L)

    # Nearest neighbors
    NN1 <- apply(d, 1L, which.min)

    # Change order and margin for nearest neighbor search
    d2 <- dtw_lb(data_reinterpolated[6:50], data_reinterpolated[1:5],
                 window.size = 20L, nn.margin = 2L)

    # Nearest neighbors *column-wise*
    NN2 <- apply(d2, 2L, which.min)

    # Same results?
    expect_identical(NN1, NN2, info = "Indices of nearest neighbors")
})

# ==================================================================================================
# dtw_basic
# ==================================================================================================

test_that("Backtracking in dtw_basic() works.", {
    univariate_result <- dtw_basic(x_uv, x_uv, backtrack = TRUE)
    expect_identical(univariate_result$index1, univariate_result$index2)
    multivariate_result <- dtw_basic(x_mv, x_mv, backtrack = TRUE)
    expect_identical(multivariate_result$index1, multivariate_result$index2)
})

# ==================================================================================================
# GAK
# ==================================================================================================

test_that("GAK can estimate sigma.", {
    gak_distance <- GAK(data[[1L]], data[[100L]])
    expect_gt(attr(gak_distance, "sigma"), 0)
})

# ==================================================================================================
# sdtw
# ==================================================================================================

test_that("Inconsistencies in parameters for sdtw() are detected.", {
    expect_error(sdtw(x_uv, x_uv, gamma = -0.01), "gamma.*positive")
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
