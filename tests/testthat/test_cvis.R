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

     cvis <- c(internal_cvis, external_cvis)

     expect_true(all(replicate(100L, {
          considered_cvis <- sample(cvis, sample(length(cvis), 1L))
          this_cvis <- cvi(pc_mv, rep(1:4, each = 5), considered_cvis)
          all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
     })))

     skip_on_cran()

     expect_equal_to_reference(base_cvis, file_name(base_cvis))
})

# =================================================================================================
# external
# =================================================================================================

test_that("external CVI calculations are consistent regardless of quantity or order of CVIs computed", {
     base_cvis <- cvi(labels_shuffled, CharTrajLabels, "external")

     expect_true(all(replicate(1000L, {
          considered_cvis <- sample(external_cvis, sample(length(external_cvis), 1L))
          this_cvis <- cvi(labels_shuffled, CharTrajLabels, considered_cvis)
          all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
     })))
})
