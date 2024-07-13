# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# multiple k
# ==================================================================================================

test_that("Multiple k works as expected.", {
    hc_k <- tsclust(data_reinterpolated, type = "h", k = 2L:5L, distance = "L2", seed = 938)
    expect_identical(length(hc_k), 4L)
    hc_k <- lapply(hc_k, reset_nondeterministic)
    assign("hc_k", hc_k, persistent)
})

# ==================================================================================================
# hierarchical algorithms
# ==================================================================================================

test_that("Hierarchical clustering works as expected.", {
    ## ---------------------------------------------------------- all
    hc_all <- tsclust(data, type = "hierarchical", k = 20L,
                      distance = "sbd",
                      control = hierarchical_control(method = "all"))
    hc_all <- lapply(hc_all, reset_nondeterministic)
    assign("hc_all", hc_all, persistent)

    ## ---------------------------------------------------------- with provided distmat
    id_avg <- which(sapply(hc_all, slot, "method") == "average")
    distmat <- as.matrix(hc_all[[1L]]@distmat)
    attr(distmat, "method") <- attr(hc_all[[1L]]@distmat, "method")
    expect_output(
        hc_avg <- tsclust(data, type = "hierarchical", k = 20L,
                          distance = "sbd", trace = TRUE,
                          control = hierarchical_control(method = "average",
                                                         distmat = distmat)),
        "provided"
    )
    expect_identical(hc_all[[id_avg]]@cluster, hc_avg@cluster)
    expect_identical(hc_all[[id_avg]]@centroids, hc_avg@centroids)
    expect_identical(hc_avg@distance, "SBD")

    ## ---------------------------------------------------------- errors with provided distmat
    expect_error(
        hc_avg <- tsclust(data, type = "hierarchical", k = 20L,
                          distance = "sbd",
                          control = hierarchical_control(method = "average",
                                                         distmat = distmat[1L:2L, 1L:2L])),
        "distance matrix"
    )

    attr(distmat, "method") <- NULL
    expect_error(
        hc_avg <- tsclust(data, type = "hierarchical", k = 20L,
                          distance = "sbd",
                          control = hierarchical_control(method = "average",
                                                         distmat = distmat)),
        "'method' attribute"
    )

    ## ---------------------------------------------------------- non-symmetric
    expect_warning({
        hc_lbi <- tsclust(data_reinterpolated, type = "hierarchical", k = 20L,
                          distance = "lbi",
                          control = hierarchical_control(method = "all"),
                          args = tsclust_args(dist = list(window.size = 17L)))
    })
    hc_lbi <- lapply(hc_lbi, reset_nondeterministic)
    assign("hc_lbi", hc_lbi, persistent)

    ## ---------------------------------------------------------- custom centroid
    hc_cent <- tsclust(data, type = "hierarchical", k = 20L,
                       distance = "sbd",
                       preproc = zscore, centroid = shape_extraction,
                       seed = 320,
                       control = hierarchical_control(method = "all"))
    hc_cent <- lapply(hc_cent, reset_nondeterministic)
    assign("hc_cent", hc_cent, persistent)

    hc_cent2 <- tsclust(data_subset, type = "hierarchical", k = 2L,
                        distance = "sbd", centroid = sdtw_cent,
                        seed = 320)
    hc_cent2 <- reset_nondeterministic(hc_cent2)
    expect_identical(hc_cent2@centroid, "sdtw_cent")
    assign("hc_cent2", hc_cent2, persistent)
})

# ==================================================================================================
# cumstom hierarchical function
# ==================================================================================================

test_that("A valid custom hierarchical function works as expected.", {
    suppressPackageStartupMessages(require(cluster))

    hc_diana <- tsclust(data, type = "hierarchical", k = 20L,
                        distance = "sbd",
                        control = hierarchical_control(method = diana))
    expect_s4_class(hc_diana, "HierarchicalTSClusters")
    hc_diana <- reset_nondeterministic(hc_diana)
    hc_diana$call <- NULL
    assign("hc_diana", hc_diana, persistent)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
