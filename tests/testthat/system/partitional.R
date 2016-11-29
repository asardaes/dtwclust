context("\tPartitional")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# multiple k and repetitions
# =================================================================================================

test_that("Multiple k and multiple repetitions work as expected.", {
    pc_k <- dtwclust(data_reinterpolated, type = "p", k = 2L:5L,
                     distance = "L2", centroid = "pam",
                     seed = 938)

    expect_identical(length(pc_k), 4L)

    pc_rep <- dtwclust(data_reinterpolated, type = "p", k = 20L,
                       distance = "L2", centroid = "pam",
                       seed = 938, control = list(nrep = 2L))

    expect_identical(length(pc_rep), 2L)

    pc_krep <- dtwclust(data_reinterpolated, type = "p", k = 20L:22L,
                        distance = "L2", centroid = "pam",
                        seed = 938, control = list(nrep = 2L))

    expect_identical(length(pc_krep), 6L)

    pc_k <- lapply(pc_k, reset_nondeterministic)
    pc_rep <- lapply(pc_rep, reset_nondeterministic)
    pc_krep <- lapply(pc_krep, reset_nondeterministic)

    pc_rep <- lapply(pc_rep, function(obj) {
        obj@call <- as.call(list("zas", a = 1))
    })

    pc_krep <- lapply(pc_krep, function(obj) {
        obj@call <- as.call(list("zas", a = 1))
    })

    expect_identical(pc_rep[1L:2L], pc_krep[1L:2L])

    ## refs
    assign("pc_k", pc_k, persistent)
    assign("pc_rep", pc_rep, persistent)
    assign("pc_krep", pc_krep, persistent)
})

# =================================================================================================
# partitional algorithms
# =================================================================================================

test_that("Partitional clustering works as expected.", {
    ## ---------------------------------------------------------- dtw
    pc_dtwb <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                        distance = "dtw_basic", centroid = "pam",
                        seed = 938, control = list(window.size = 20L))

    pc_dtwb_npampre <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam",
                                seed = 938, control = list(window.size = 20L,
                                                           pam.precompute = FALSE))

    pc_dtwb_distmat <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam",
                                seed = 938, control = list(window.size = 20L),
                                distmat = pc_dtwb@distmat)

    pc_dtwlb <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                         distance = "dtw_lb", centroid = "pam",
                         seed = 938, control = list(window.size = 20L,
                                                    pam.precompute = FALSE))

    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre@cluster, info = "pam.pre vs no pam.pre")
    expect_identical(pc_dtwb@cluster, pc_dtwb_distmat@cluster, info = "distmat supplied")
    expect_identical(pc_dtwb@cluster, pc_dtwlb@cluster, info = "dtw_basic vs dtw_lb")

    pc_dtwb <- reset_nondeterministic(pc_dtwb)
    pc_dtwb_npampre <- reset_nondeterministic(pc_dtwb_npampre)
    pc_dtwb_distmat <- reset_nondeterministic(pc_dtwb_distmat)
    pc_dtwlb <- reset_nondeterministic(pc_dtwlb)

    assign("pc_dtwb", pc_dtwb, persistent)
    assign("pc_dtwb_npampre", pc_dtwb_npampre, persistent)
    assign("pc_dtwb_distmat", pc_dtwb_distmat, persistent)
    assign("pc_dtwlb", pc_dtwlb, persistent)

    ## ---------------------------------------------------------- k-Shape
    pc_kshape <- dtwclust(data_subset, type = "p", k = 4L,
                          distance = "sbd", centroid = "shape",
                          seed = 938)

    pc_kshape <- reset_nondeterministic(pc_kshape)

    assign("pc_kshape", pc_kshape, persistent)

    ## ---------------------------------------------------------- dba
    pc_dba <- dtwclust(data_subset, type = "p", k = 4L,
                       distance = "dtw_basic", centroid = "dba",
                       seed = 938)

    pc_dba <- reset_nondeterministic(pc_dba)

    assign("pc_dba", pc_dba, persistent)

    ## ---------------------------------------------------------- multivariate pam
    pc_mv_pam <- dtwclust(data_multivariate, type = "p", k = 4L,
                          distance = "dtw_basic", centroid = "pam",
                          seed = 938)

    pc_mv_pam <- reset_nondeterministic(pc_mv_pam)

    assign("pc_mv_pam", pc_mv_pam, persistent)

    ## ---------------------------------------------------------- multivariate dba
    pc_mv_dba <- dtwclust(data_multivariate, type = "p", k = 4L,
                          distance = "dtw_basic", centroid = "dba",
                          seed = 938)

    pc_mv_dba <- reset_nondeterministic(pc_mv_dba)

    assign("pc_mv_dba", pc_mv_dba, persistent)
})

# =================================================================================================
# TADPole
# =================================================================================================

test_that("TADPole works as expected", {
    ## ---------------------------------------------------------- TADPole
    pc_tadp <- dtwclust(data_reinterpolated_subset, type = "t", k = 4L,
                        dc = 1.5, control = list(window.size = 20L),
                        seed = 938)

    pc_tadp <- reset_nondeterministic(pc_tadp)

    assign("pc_tadp", pc_tadp, persistent)

    ## ---------------------------------------------------------- TADPole with LBI
    pc_tadp_lbi <- dtwclust(data_reinterpolated_subset, type = "t", k = 4L,
                            dc = 1.5, control = list(window.size = 20L),
                            lb = "lbi", seed = 938)

    pc_tadp_lbi <- reset_nondeterministic(pc_tadp_lbi)

    assign("pc_tadp_lbi", pc_tadp_lbi, persistent)

    ## ---------------------------------------------------------- TADPole with custom centroid
    pc_tadp_cent <- dtwclust(data_reinterpolated_subset, type = "t", k = 4L,
                             preproc = zscore, centroid = shape_extraction,
                             dc = 1.5, control = list(window.size = 20L),
                             seed = 938)

    pc_tadp_cent <- reset_nondeterministic(pc_tadp_cent)

    assign("pc_tadp_cent", pc_tadp_cent, persistent)
})

# =================================================================================================
# cluster reinitialization
# =================================================================================================

## this case causes clusters to become empty
test_that("Cluster reinitialization in partitional dtwclust works.", {
    suppressWarnings(pc_cr <- dtwclust(data_reinterpolated, k = 20,
                                       distance = "lbk", centroid = "mean",
                                       seed=31231,
                                       control = list(window.size = 19L,
                                                      iter.max = 10L)))

    expect_false(pc_cr@converged)

    pc_cr <- reset_nondeterministic(pc_cr)
    assign("pc_cr", pc_cr, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
