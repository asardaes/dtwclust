context("Test only CVIs")

# =================================================================================================
# both
# =================================================================================================

pc_mv <- dtwclust(data_multivariate, type = "partitional", k = 4,
                  distance = "dtw", centroid = "pam",
                  preproc = NULL, control = list(window.size = 18L), seed = 123,
                  dist.method = "L1")

base_cvis <- cvi(pc_mv, rep(1:4, each = 5), "valid")

cvis <- c(internal_cvis, external_cvis)

test_that("dtwclust CVI calculations are consistent regardless of quantity or order of CVIs computed",
          expect_true(all(replicate(100L, {
               considered_cvis <- sample(cvis, sample(length(cvis), 1L))
               this_cvis <- cvi(pc_mv, rep(1:4, each = 5), considered_cvis)
               all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
          }))))

# =================================================================================================
# external
# =================================================================================================

base_cvis <- cvi(labels_shuffled, CharTrajLabels, "external")

test_that("external CVI calculations are consistent regardless of quantity or order of CVIs computed",
          expect_true(all(replicate(1000L, {
               considered_cvis <- sample(external_cvis, sample(length(external_cvis), 1L))
               this_cvis <- cvi(labels_shuffled, CharTrajLabels, considered_cvis)
               all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
          }))))
