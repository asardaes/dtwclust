context("\tHierarchical")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# multiple k
# =================================================================================================

test_that("Multiple k works as expected.", {
    hc_k <- dtwclust(data_reinterpolated, type = "h", k = 2L:5L,
                     distance = "L2", seed = 938)

    hc_k2 <- tsclust(data_reinterpolated, type = "h", k = 2L:5L,
                     distance = "L2", seed = 938)

    for (i in seq_along(hc_k)) {
        expect_identical(hc_k[[i]]@cluster, hc_k2[[i]]@cluster)
        expect_identical(hc_k[[i]]@centroids, hc_k2[[i]]@centroids)
        expect_identical(hc_k[[i]]@distmat, hc_k2[[i]]@distmat)
    }

    expect_identical(length(hc_k), 4L)

    skip_on_cran()

    hc_k <- lapply(hc_k, reset_nondeterministic)
    assign("hc_k", hc_k, persistent)
})

# =================================================================================================
# hierarchical algorithms
# =================================================================================================

test_that("Hierarchical clustering works as expected.", {
    ## ---------------------------------------------------------- all
    hc_all <- dtwclust(data, type = "hierarchical", k = 20L,
                       distance = "sbd", method = "all")

    hc_all2 <- tsclust(data, type = "hierarchical", k = 20L,
                       distance = "sbd",
                       control = hierarchical_control(method = "all"))

    for (i in seq_along(hc_all)) {
        expect_identical(hc_all[[i]]@cluster, hc_all2[[i]]@cluster)
        expect_identical(hc_all[[i]]@centroids, hc_all2[[i]]@centroids)
        expect_identical(hc_all[[i]]@distmat, hc_all2[[i]]@distmat)
    }

    hc_all <- lapply(hc_all, reset_nondeterministic)

    assign("hc_all", hc_all, persistent)

    ## ---------------------------------------------------------- non-symmetric
    hc_lbi <- dtwclust(data_reinterpolated, type = "hierarchical", k = 20L,
                       distance = "lbi", method = "all",
                       control = list(window.size = 17L))

    hc_lbi2 <- tsclust(data_reinterpolated, type = "hierarchical", k = 20L,
                       distance = "lbi",
                       control = hierarchical_control(method = "all"),
                       args = tsclust_args(dist = list(window.size = 17L)))

    for (i in seq_along(hc_lbi)) {
        expect_identical(hc_lbi[[i]]@cluster, hc_lbi2[[i]]@cluster)
        expect_identical(hc_lbi[[i]]@centroids, hc_lbi2[[i]]@centroids)
        expect_identical(hc_lbi[[i]]@distmat, hc_lbi2[[i]]@distmat)
    }

    hc_lbi <- lapply(hc_lbi, reset_nondeterministic)

    assign("hc_lbi", hc_lbi, persistent)

    ## ---------------------------------------------------------- custom centroid
    hc_cent <- dtwclust(data, type = "hierarchical", k = 20L,
                        distance = "sbd", method = "all",
                        preproc = zscore, centroid = shape_extraction,
                        seed = 320)

    hc_cent2 <- tsclust(data, type = "hierarchical", k = 20L,
                        distance = "sbd",
                        preproc = zscore, centroid = shape_extraction,
                        seed = 320,
                        control = hierarchical_control(method = "all"))

    for (i in seq_along(hc_cent)) {
        expect_identical(hc_cent[[i]]@cluster, hc_cent2[[i]]@cluster)
        expect_identical(hc_cent[[i]]@centroids, hc_cent2[[i]]@centroids)
        expect_identical(hc_cent[[i]]@distmat, hc_cent2[[i]]@distmat)
    }

    hc_cent <- lapply(hc_cent, reset_nondeterministic)

    assign("hc_cent", hc_cent, persistent)
})

# =================================================================================================
# cumstom hierarchical function
# =================================================================================================

test_that("A valid custom hierarchical function works as expected.", {
    require(cluster)

    hc_diana <- dtwclust(data, type = "hierarchical", k = 20L,
                         distance = "sbd", method = diana)

    hc_diana2 <- tsclust(data, type = "hierarchical", k = 20L,
                         distance = "sbd",
                         control = hierarchical_control(method = diana))

    expect_identical(hc_diana@cluster, hc_diana2@cluster)
    expect_identical(hc_diana@centroids, hc_diana2@centroids)
    expect_identical(hc_diana@distmat, hc_diana2@distmat)

    hc_diana <- reset_nondeterministic(hc_diana)
    hc_diana$call <- NULL

    assign("hc_diana", hc_diana, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
