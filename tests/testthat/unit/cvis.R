context("\tCVIs")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

internal_cvis <- c("Sil", "D", "DB", "DBstar", "CH", "SF", "COP")
external_cvis <- c("RI", "ARI", "J", "FM", "VI")

# =================================================================================================
# both internal and external
# =================================================================================================

test_that("CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    pc_mv <- dtwclust(data_multivariate, type = "partitional", k = 4L,
                      distance = "dtw_basic", centroid = "pam",
                      preproc = NULL, control = list(window.size = 18L), seed = 123,
                      dist.method = "L1")

    pc_mv2 <- as(pc_mv, "TSClusters")

    base_cvis <- cvi(pc_mv, rep(1L:4L, each = 5L), "valid")
    i_cvis <- cvi(pc_mv, type = "internal")
    e_cvis <- cvi(pc_mv, rep(1L:4L, each = 5L), type = "external")
    expect_identical(base_cvis, c(e_cvis, i_cvis))

    cvis <- c(internal_cvis, external_cvis)

    expect_true(all(
        times(50L) %dopar% {
            considered_cvis <- sample(cvis, sample(length(cvis), 1L))
            this_cvis <- cvi(pc_mv, rep(1L:4L, each = 5L), considered_cvis)
            this_cvis2 <- cvi(pc_mv2, rep(1L:4L, each = 5L), considered_cvis)
            all(base_cvis[considered_cvis] == this_cvis[considered_cvis]) &&
                identical(this_cvis, this_cvis2)
        }
    ),
    info = "A random number of CVIs are calculated and compared against the base ones, and should always be equal.")

    ## when missing elements
    pc_mv2@distmat <- pc_mv@distmat <- NULL
    this_cvis <- cvi(pc_mv, type = "internal")
    this_cvis2 <- cvi(pc_mv2, type = "internal")
    considered_cvis <- names(this_cvis)
    expect_true(all(base_cvis[considered_cvis] == this_cvis))
    expect_identical(this_cvis, this_cvis2)

    pc_mv2@datalist <- pc_mv@datalist <- list()
    expect_warning(this_cvis <- cvi(pc_mv, type = "internal"))
    expect_warning(this_cvis2 <- cvi(pc_mv2, type = "internal"))
    considered_cvis <- names(this_cvis)
    expect_true(all(base_cvis[considered_cvis] == this_cvis))
    expect_identical(this_cvis, this_cvis2)

    ## refs
    assign("base_cvis", base_cvis, persistent)
})

# =================================================================================================
# external
# =================================================================================================

test_that("external CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    expect_error(cvi(labels_shuffled, type = "external"))
    base_cvis <- cvi(labels_shuffled, CharTrajLabels, "external")
    export <- c("external_cvis", "labels_shuffled", "CharTrajLabels")

    expect_true(all(
        foreach(iterators::icount(1000L), .combine = c, .export = export) %dopar% {
            considered_cvis <- sample(external_cvis, sample(length(external_cvis), 1L))
            this_cvis <- cvi(labels_shuffled, CharTrajLabels, considered_cvis)
            all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
        }
    ),
    info = "A random number of CVIs are calculated and compared against the base ones, and should always be equal.")
})

# =================================================================================================
# hierarchical/tadpole/fuzzy cases
# =================================================================================================

test_that("CVIs work also for hierarchical and TADPole", {
    tadp <- dtwclust(data_reinterpolated_subset, type = "t", k = 4L,
                     dc = 1.5, control = list(window.size = 18L))

    hc <- dtwclust(data_reinterpolated_subset, type = "h", k = 4L,
                   distance = "gak", sigma = 100,
                   control = list(window.size = 18L))

    fc <- dtwclust(data_reinterpolated_subset, type = "f", k = 4L, distance = "L2")

    expect_error(cvi(fc, labels_subset))
    cvis_tadp <- cvi(tadp, labels_subset)
    cvis_hc <- cvi(hc, labels_subset)
    cvis_tadp2 <- cvi(as(tadp, "TSClusters"), labels_subset)
    cvis_hc2 <- cvi(as(hc, "TSClusters"), labels_subset)

    expect_error(cvi(as(fc, "TSClusters"), labels_subset))
    expect_identical(cvis_tadp2, cvis_tadp)
    expect_identical(cvis_hc2, cvis_hc)

    ## refs
    assign("cvis_tadp", cvis_tadp, persistent)
    assign("cvis_hc", cvis_hc, persistent)
})

test_that("CVIs work also for hierarchical and TADPole with custom centroid", {
    tadp <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                    centroid = shape_extraction,
                    control = tadpole_control(dc = 1.5, window.size = 18L))

    hc <- tsclust(data_reinterpolated_subset, type = "h", k = 4L,
                  distance = "sbd", centroid = shape_extraction)

    cvis_tadp_cent <- cvi(tadp, labels_subset)
    cvis_hc_cent <- cvi(hc, labels_subset)

    ## refs
    assign("cvis_tadp_cent", cvis_tadp_cent, persistent)
    assign("cvis_hc_cent", cvis_hc_cent, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
