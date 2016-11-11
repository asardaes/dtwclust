context("Test CVIs")

# =================================================================================================
# both internal and external
# =================================================================================================

test_that("dtwclust CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    pc_mv <- dtwclust(data_multivariate, type = "partitional", k = 4,
                      distance = "dtw_basic", centroid = "pam",
                      preproc = NULL, control = list(window.size = 18L), seed = 123,
                      dist.method = "L1")

    base_cvis <- cvi(pc_mv, rep(1:4, each = 5), "valid")
    i_cvis <- cvi(pc_mv, type = "internal")
    e_cvis <- cvi(pc_mv, rep(1:4, each = 5), type = "external")
    expect_identical(base_cvis, c(e_cvis, i_cvis))

    cvis <- c(internal_cvis, external_cvis)

    expect_true(all(replicate(100L, {
        considered_cvis <- sample(cvis, sample(length(cvis), 1L))
        this_cvis <- cvi(pc_mv, rep(1:4, each = 5), considered_cvis)
        all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
    })))

    pc_mv@datalist <- list()
    pc_mv@distmat <- NULL
    expect_warning(cvi(pc_mv, type = "valid"))

    skip_on_cran()

    expect_equal_to_reference(base_cvis, file_name(base_cvis))
})

# =================================================================================================
# external
# =================================================================================================

test_that("external CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    expect_error(cvi(labels_shuffled, type = "external"))

    base_cvis <- cvi(labels_shuffled, CharTrajLabels, "external")

    expect_true(all(replicate(1000L, {
        considered_cvis <- sample(external_cvis, sample(length(external_cvis), 1L))
        this_cvis <- cvi(labels_shuffled, CharTrajLabels, considered_cvis)
        all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
    })))
})
