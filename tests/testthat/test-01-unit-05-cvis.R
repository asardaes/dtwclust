# ==================================================================================================
# setup
# ==================================================================================================

# Original objects in env
ols <- ls()

internal_cvis <- c("Sil", "D", "DB", "DBstar", "CH", "SF", "COP")
external_cvis <- c("RI", "ARI", "J", "FM", "VI")
internal_fuzzy_cvis <- c("MPC", "K", "T", "SC", "PBMF")
external_fuzzy_cvis <- c("ARI", "RI", "VI", "NMIM")

# ==================================================================================================
# both internal and external
# ==================================================================================================

test_that("CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    pc_mv <- tsclust(data_multivariate, type = "partitional", k = 4L,
                     distance = "dtw_basic", centroid = "pam",
                     args = tsclust_args(dist = list(window.size = 18L)),
                     seed = 123)

    expect_warning(
        expect_warning(
            base_cvis <- cvi(pc_mv, rep(1L:4L, each = 5L), "valid")
        )
    )
    expect_warning(
        expect_warning(
            i_cvis <- cvi(pc_mv, type = "internal")
        )
    )
    e_cvis <- cvi(pc_mv, rep(1L:4L, each = 5L), type = "external")
    expect_identical(base_cvis, c(e_cvis, i_cvis))

    cvis <- c(internal_cvis, external_cvis)

    expect_true(all(
        times(50L) %dopar% {
            considered_cvis <- sample(cvis, sample(length(cvis), 1L))
            suppressWarnings(
                this_cvis <- cvi(pc_mv, rep(1L:4L, each = 5L), considered_cvis)
            )
            all(base_cvis[considered_cvis] == this_cvis[considered_cvis])
        }
    ),
    info = paste("A random number of CVIs are calculated and compared against the base ones,",
                 "and should always be equal."))

    # when missing elements
    pc_mv@distmat <- NULL
    expect_warning(
        expect_warning(
            this_cvis <- cvi(pc_mv, type = "internal")
        )
    )
    considered_cvis <- names(this_cvis)
    expect_true(all(base_cvis[considered_cvis] == this_cvis))

    pc_mv@datalist <- list()
    expect_warning(
        expect_warning(
            expect_warning(
                this_cvis <- cvi(pc_mv, type = "internal")
            )
        )
    )
    considered_cvis <- names(this_cvis)
    expect_true(all(base_cvis[considered_cvis] == this_cvis))

    # refs
    assign("base_cvis", base_cvis, persistent)
})

# ==================================================================================================
# external
# ==================================================================================================

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
    info = paste("A random number of CVIs are calculated and compared against the base ones,",
                 "and should always be equal."))
})

# ==================================================================================================
# fuzzy
# ==================================================================================================

fc <- tsclust(data_subset, "f", 4L, distance = "sbd", centroid = "fcmdd", seed = 32890L)
base_fcvis <- cvi(fc, labels_subset, type = "valid")

test_that("Internal and external fuzzy CVIs are identical to the valid ones", {
    icvis <- cvi(fc, type = "internal")
    ecvis <- cvi(fc, labels_subset, type = "external")
    expect_identical(icvis[internal_fuzzy_cvis], base_fcvis[internal_fuzzy_cvis])
    expect_identical(ecvis[external_fuzzy_cvis], base_fcvis[external_fuzzy_cvis])
    expect_identical(base_fcvis[external_fuzzy_cvis], cvi(fc@fcluster, labels_subset))
    expect_error(cvi(fc@fcluster, type = "internal"))
    expect_error(cvi(fc, type = "external"))
})

# Internal
test_that("Internal fuzzy CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    # parallel times() below won't detect 'fc' otherwise -.-
    fc <- fc

    internal_fcvis <- base_fcvis[internal_fuzzy_cvis]
    cvis <- internal_fuzzy_cvis
    `%op%` <- dtwclust:::`%op%` # avoid stupid parallel warnings

    expect_true(all(
        times(50L) %op% {
            considered_cvis <- sample(cvis, sample(length(cvis), 1L))
            this_cvis <- cvi(fc, type = considered_cvis)
            all(internal_fcvis[considered_cvis] == this_cvis[considered_cvis])
        }
    ),
    info = paste("A random number of internal fuzzy CVIs are calculated and compared against the base ones,",
                 "and should always be equal."))

    # when missing elements
    fc@datalist <- list()
    expect_warning(this_cvis <- cvi(fc))
    considered_cvis <- names(this_cvis)
    expect_true(all(internal_fcvis[considered_cvis] == this_cvis))

    # refs
    assign("internal_fcvis", internal_fcvis, persistent)
})

# External
test_that("External fuzzy CVI calculations are consistent regardless of quantity or order of CVIs computed", {
    # parallel times() below won't detect 'fc' otherwise -.-
    fc <- fc

    external_fcvis <- base_fcvis[external_fuzzy_cvis]
    cvis <- external_fuzzy_cvis
    `%op%` <- dtwclust:::`%op%` # avoid stupid parallel warnings

    expect_true(all(
        times(50L) %op% {
            considered_cvis <- sample(cvis, sample(length(cvis), 1L))
            this_cvis <- cvi(fc, rep(1L:4L, each = 5L), type = considered_cvis)
            all(external_fcvis[considered_cvis] == this_cvis[considered_cvis])
        }
    ),
    info = paste("A random number of external fuzzy CVIs are calculated and compared against the base ones,",
                 "and should always be equal."))

    # refs
    assign("external_fcvis", external_fcvis, persistent)
})

# ==================================================================================================
# hierarchical/tadpole cases
# ==================================================================================================

test_that("CVIs work also for hierarchical and TADPole", {
    tadp <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                    control = tadpole_control(1.5, 18L))

    hc <- tsclust(data_reinterpolated_subset, type = "h", k = 4L,
                  distance = "gak", sigma = 100,
                  window.size = 18L)

    expect_warning(
        expect_warning(
            cvis_tadp <- cvi(tadp, labels_subset)
        )
    )
    cvis_hc <- cvi(hc, labels_subset)

    # refs
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
    expect_length(cvis_tadp_cent, 12L)
    cvis_hc_cent <- cvi(hc, labels_subset)
    expect_length(cvis_hc_cent, 12L)

    # refs
    assign("cvis_tadp_cent", cvis_tadp_cent, persistent)
    assign("cvis_hc_cent", cvis_hc_cent, persistent)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
