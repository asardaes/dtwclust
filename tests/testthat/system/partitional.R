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

    pc_k2 <- tsclust(data_reinterpolated, type = "p", k = 2L:5L,
                     distance = "L2", centroid = "pam",
                     seed = 938)

    for (i in seq_along(pc_k)) {
        expect_identical(pc_k[[i]]@cluster, pc_k2[[i]]@cluster)
        expect_identical(pc_k[[i]]@centroids, pc_k2[[i]]@centroids)
        expect_identical(pc_k[[i]]@distmat, pc_k2[[i]]@distmat)
    }

    pc_rep <- dtwclust(data_reinterpolated, type = "p", k = 20L,
                       distance = "L2", centroid = "pam",
                       seed = 938, control = list(nrep = 2L))

    expect_identical(length(pc_rep), 2L)

    pc_rep2 <- tsclust(data_reinterpolated, type = "p", k = 20L,
                       distance = "L2", centroid = "pam",
                       seed = 938, control = partitional_control(nrep = 2L))

    for (i in seq_along(pc_rep)) {
        expect_identical(pc_rep[[i]]@cluster, pc_rep2[[i]]@cluster)
        expect_identical(pc_rep[[i]]@centroids, pc_rep2[[i]]@centroids)
        expect_identical(pc_rep[[i]]@distmat, pc_rep2[[i]]@distmat)
    }

    pc_krep <- dtwclust(data_reinterpolated, type = "p", k = 20L:22L,
                        distance = "L2", centroid = "pam",
                        seed = 938, control = list(nrep = 2L))

    expect_identical(length(pc_krep), 6L)

    pc_krep2 <- tsclust(data_reinterpolated, type = "p", k = 20L:22L,
                        distance = "L2", centroid = "pam",
                        seed = 938, control = partitional_control(nrep = 2L))

    for (i in seq_along(pc_krep)) {
        expect_identical(pc_krep[[i]]@cluster, pc_krep2[[i]]@cluster)
        expect_identical(pc_krep[[i]]@centroids, pc_krep2[[i]]@centroids)
        expect_identical(pc_krep[[i]]@distmat, pc_krep2[[i]]@distmat)
    }

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

    pc_dtwb2 <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                        distance = "dtw_basic", centroid = "pam", seed = 938,
                        args = tsclust_args(dist = list(window.size = 20L)))

    pc_dtwb_npampre <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam",
                                seed = 938, control = list(window.size = 20L,
                                                           pam.precompute = FALSE))

    pc_dtwb_npampre2 <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam", seed = 938,
                                args = tsclust_args(dist = list(window.size = 20L)),
                                control = partitional_control(pam.precompute = FALSE))

    pc_dtwb_npampre3 <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam", seed = 938,
                                args = tsclust_args(dist = list(window.size = 20L)),
                                control = partitional_control(pam.precompute = FALSE,
                                                              pam.sparse = FALSE))

    pc_dtwb_distmat <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam",
                                seed = 938, control = list(window.size = 20L),
                                distmat = pc_dtwb@distmat)

    pc_dtwb_distmat2 <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam", seed = 938,
                                args = tsclust_args(dist = list(window.size = 20L)),
                                control = partitional_control(distmat = pc_dtwb@distmat))

    pc_dtwlb <- dtwclust(data_reinterpolated_subset, type = "p", k = 4L,
                         distance = "dtw_lb", centroid = "pam",
                         seed = 938, control = list(window.size = 20L,
                                                    pam.precompute = FALSE))

    pc_dtwlb2 <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                         distance = "dtw_lb", centroid = "pam", seed = 938,
                         args = tsclust_args(dist = list(window.size = 20L)),
                         control = partitional_control(pam.precompute = FALSE))

    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre@cluster, info = "pam.pre vs no pam.pre")
    expect_identical(pc_dtwb@cluster, pc_dtwb_distmat@cluster, info = "distmat supplied")
    expect_identical(pc_dtwb@cluster, pc_dtwlb@cluster, info = "dtw_basic vs dtw_lb")
    expect_identical(pc_dtwb@cluster, pc_dtwb2@cluster, info = "TSC")
    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre2@cluster, info = "TSC no pam.pre")
    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre3@cluster, info = "TSC no pam.pre nor sparse")
    expect_identical(pc_dtwb@cluster, pc_dtwb_distmat2@cluster, info = "TSC distmat supplied")
    expect_identical(pc_dtwb@cluster, pc_dtwlb2@cluster, info = "TSC dtw_basic vs dtw_lb")

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

    pc_kshape2 <- tsclust(data_subset, type = "p", k = 4L,
                          distance = "sbd", centroid = "shape",
                          seed = 938)

    expect_identical(pc_kshape@cluster, pc_kshape2@cluster)
    expect_identical(pc_kshape@centroids, pc_kshape2@centroids)

    pc_kshape <- reset_nondeterministic(pc_kshape)

    assign("pc_kshape", pc_kshape, persistent)

    ## ---------------------------------------------------------- dba
    pc_dba <- dtwclust(data_subset, type = "p", k = 4L,
                       distance = "dtw_basic", centroid = "dba",
                       seed = 938)

    pc_dba2 <- tsclust(data_subset, type = "p", k = 4L,
                       distance = "dtw_basic", centroid = "dba",
                       args = tsclust_args(cent = list(max.iter = 15L)),
                       seed = 938)

    expect_identical(pc_dba@cluster, pc_dba2@cluster)
    expect_identical(pc_dba@centroids, pc_dba2@centroids)

    pc_dba <- reset_nondeterministic(pc_dba)

    assign("pc_dba", pc_dba, persistent)

    ## ---------------------------------------------------------- multivariate pam
    pc_mv_pam <- dtwclust(data_multivariate, type = "p", k = 4L,
                          distance = "dtw_basic", centroid = "pam",
                          seed = 938)

    pc_mv_pam2 <- tsclust(data_multivariate, type = "p", k = 4L,
                          distance = "dtw_basic", centroid = "pam",
                          seed = 938)

    expect_identical(pc_mv_pam@cluster, pc_mv_pam2@cluster)
    expect_identical(pc_mv_pam@centroids, pc_mv_pam2@centroids)

    pc_mv_pam <- reset_nondeterministic(pc_mv_pam)

    assign("pc_mv_pam", pc_mv_pam, persistent)

    ## ---------------------------------------------------------- multivariate dba
    pc_mv_dba <- dtwclust(data_multivariate, type = "p", k = 4L,
                          distance = "dtw_basic", centroid = "dba",
                          seed = 938)

    pc_mv_dba2 <- tsclust(data_multivariate, type = "p", k = 4L,
                          distance = "dtw_basic", centroid = "dba",
                          args = tsclust_args(cent = list(max.iter = 15L)),
                          seed = 938)

    expect_identical(pc_mv_dba@cluster, pc_mv_dba2@cluster)
    expect_identical(pc_mv_dba@centroids, pc_mv_dba2@centroids)

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

    pc_tadp2 <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                        control = tadpole_control(dc = 1.5, window.size = 20L),
                        seed = 938)

    expect_identical(pc_tadp@cluster, pc_tadp2@cluster)
    expect_identical(pc_tadp@centroids, pc_tadp2@centroids)

    pc_tadp <- reset_nondeterministic(pc_tadp)

    assign("pc_tadp", pc_tadp, persistent)

    ## ---------------------------------------------------------- TADPole with LBI
    pc_tadp_lbi <- dtwclust(data_reinterpolated_subset, type = "t", k = 4L,
                            dc = 1.5, control = list(window.size = 20L),
                            lb = "lbi", seed = 938)

    pc_tadp_lbi2 <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                            control = tadpole_control(dc = 1.5,
                                                      window.size = 20L,
                                                      lb = "lbi"),
                            seed = 938)

    expect_identical(pc_tadp_lbi@cluster, pc_tadp_lbi2@cluster)
    expect_identical(pc_tadp_lbi@centroids, pc_tadp_lbi2@centroids)

    pc_tadp_lbi <- reset_nondeterministic(pc_tadp_lbi)

    assign("pc_tadp_lbi", pc_tadp_lbi, persistent)

    ## ---------------------------------------------------------- TADPole with custom centroid
    pc_tadp_cent <- dtwclust(data_reinterpolated_subset, type = "t", k = 4L,
                             preproc = zscore, centroid = shape_extraction,
                             dc = 1.5, control = list(window.size = 20L),
                             seed = 938)

    pc_tadp_cent2 <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                             preproc = zscore, centroid = shape_extraction,
                             control = tadpole_control(dc = 1.5, window.size = 20L),
                             seed = 938)

    expect_identical(pc_tadp_cent@cluster, pc_tadp_cent2@cluster)
    expect_identical(pc_tadp_cent@centroids, pc_tadp_cent2@centroids)

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
                                       seed = 31231,
                                       control = list(window.size = 19L,
                                                      iter.max = 10L)))

    expect_false(pc_cr@converged)

    suppressWarnings(pc_cr2 <- tsclust(data_reinterpolated, k = 20,
                                       distance = "lbk", centroid = "mean",
                                       seed = 31231,
                                       control = partitional_control(iter.max = 10L),
                                       args = tsclust_args(dist = list(window.size = 19L))))

    expect_false(pc_cr2@converged)

    expect_identical(pc_cr@cluster, pc_cr2@cluster)
    expect_identical(pc_cr@centroids, pc_cr2@centroids)

    pc_cr <- reset_nondeterministic(pc_cr)
    assign("pc_cr", pc_cr, persistent)

    ## test PAM too
    suppressWarnings(pc_cr <- dtwclust(data_reinterpolated, k = 20,
                                       distance = "lbk", centroid = "pam",
                                       seed = 31231,
                                       control = list(window.size = 19L,
                                                      iter.max = 10L)))

    expect_false(pc_cr@converged)

    suppressWarnings(pc_cr2 <- tsclust(data_reinterpolated, k = 20,
                                       distance = "lbk", centroid = "pam",
                                       seed = 31231,
                                       control = partitional_control(iter.max = 10L),
                                       args = tsclust_args(dist = list(window.size = 19L))))

    expect_false(pc_cr2@converged)

    expect_identical(pc_cr@cluster, pc_cr2@cluster)
    expect_identical(pc_cr@centroids, pc_cr2@centroids)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
