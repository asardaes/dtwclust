context("    Miscellaneous functions")

# ==================================================================================================
# setup
# ==================================================================================================

# Original objects in env
ols <- ls()

# ==================================================================================================
# reinterpolate
# ==================================================================================================

test_that("reinterpolate function works correctly for supported inputs.", {
    expect_length(reinterpolate(data[[1L]], 200L), 200L)
    expect_length(reinterpolate(data[[1L]], 200L, TRUE), 200L)

    expect_identical(ncol(reinterpolate(data_matrix, 200L, multivariate = FALSE)), 200L)
    expect_identical(nrow(reinterpolate(data_matrix, 200L, multivariate = TRUE)), 200L)

    data_frame <- base::as.data.frame(data_matrix)
    expect_identical(ncol(reinterpolate(data_frame, 200L, multivariate = FALSE)), 200L)
    expect_identical(nrow(reinterpolate(data_frame, 200L, multivariate = TRUE)), 200L)
})

# ==================================================================================================
# NCCc
# ==================================================================================================

test_that("NCCc function works correctly for supported inputs.", {
    x_uv <- data[[1L]]
    x_mv <- data_multivariate[[1L]]
    invalid_inputs <- list(
        "NA" = NA,
        "NULL" = NULL,
        "char" = "1",
        "bool" = TRUE,
        "empty" = numeric(),
        "miss" = c(1, NA)
    )

    expect_error(NCCc(x_uv, x_mv), regexp = "dimension")
    expect_error(NCCc(x_mv, x_uv), regexp = "dimension")
    expect_error(NCCc(x_mv, x_mv), regexp = "multivariate")
    for (input in names(invalid_inputs)) {
        expect_error(do.call("NCCc", list(x = invalid_inputs[[input]], y = x_uv)),
                     info = paste("with", input, "input in x"))
        expect_error(do.call("NCCc", list(x = x_uv, y = invalid_inputs[[input]])),
                     info = paste("with", input, "input in y"))
    }
    expect_length(NCCc(x_uv, x_uv), 2L * length(x_uv) - 1L)
})

# ==================================================================================================
# zscore
# ==================================================================================================

test_that("zscore function works correctly for supported inputs.", {
    normalized_series <- zscore(data[[1L]])
    expect_length(normalized_series, length(data[[1L]]))
    expect_true(!is.null(attributes(zscore(data[[1L]], keep.attributes = TRUE))))

    normalized_matrix <- zscore(data_matrix, multivariate = FALSE, keep.attributes = TRUE)
    expect_identical(ncol(normalized_matrix), ncol(data_matrix))
    expect_true(!is.null(attributes(normalized_matrix)))

    normalized_matrix <- zscore(data_matrix, multivariate = TRUE, keep.attributes = TRUE)
    expect_identical(nrow(normalized_matrix), nrow(data_matrix))
    expect_true(!is.null(attributes(normalized_matrix)))

    data_frame <- base::as.data.frame(data_matrix)

    normalized_data_frame <- zscore(data_frame, multivariate = FALSE, keep.attributes = TRUE)
    expect_identical(ncol(normalized_data_frame), ncol(data_frame))
    expect_true(!is.null(attributes(normalized_data_frame)))

    normalized_data_frame <- zscore(data_frame, multivariate = TRUE, keep.attributes = TRUE)
    expect_identical(nrow(normalized_data_frame), nrow(data_frame))
    expect_true(!is.null(attributes(normalized_data_frame)))
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
