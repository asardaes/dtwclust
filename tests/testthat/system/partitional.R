context("    Partitional")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# multiple k and repetitions
# ==================================================================================================

test_that("Multiple k and multiple repetitions work as expected.", {
    pc_k <- tsclust(data_reinterpolated, type = "p", k = 2L:5L,
                    distance = "L2", centroid = "pam",
                    seed = 938)
    expect_identical(length(pc_k), 4L)

    pc_rep <- tsclust(data_reinterpolated, type = "p", k = 20L,
                      distance = "L2", centroid = "pam",
                      seed = 938, control = partitional_control(nrep = 2L))
    expect_identical(length(pc_rep), 2L)

    pc_krep <- tsclust(data_reinterpolated, type = "p", k = 20L:22L,
                       distance = "L2", centroid = "pam",
                       seed = 938, control = partitional_control(nrep = 2L))
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

# ==================================================================================================
# partitional algorithms
# ==================================================================================================

test_that("Partitional clustering works as expected.", {
    ## ---------------------------------------------------------- dtw
    pc_dtwb <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                       distance = "dtw_basic", centroid = "pam", seed = 938,
                       args = tsclust_args(dist = list(window.size = 20L)))

    pc_dtwb_npampre <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                               distance = "dtw_basic", centroid = "pam", seed = 938,
                               args = tsclust_args(dist = list(window.size = 20L)),
                               control = partitional_control(pam.precompute = FALSE,
                                                             pam.sparse = TRUE))

    pc_dtwb_npampre2 <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                                distance = "dtw_basic", centroid = "pam", seed = 938,
                                args = tsclust_args(dist = list(window.size = 20L)),
                                control = partitional_control(pam.precompute = FALSE,
                                                              pam.sparse = FALSE))

    pc_dtwb_distmat <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                               distance = "dtw_basic", centroid = "pam", seed = 938,
                               args = tsclust_args(dist = list(window.size = 20L)),
                               control = partitional_control(distmat = pc_dtwb@distmat))

    pc_dtwlb <- tsclust(data_reinterpolated_subset, type = "p", k = 4L,
                        distance = "dtw_lb", centroid = "pam", seed = 938,
                        args = tsclust_args(dist = list(window.size = 20L)),
                        control = partitional_control(pam.precompute = FALSE))

    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre@cluster, info = "pam.pre vs no pam.pre")
    expect_identical(pc_dtwb@cluster, pc_dtwb_distmat@cluster, info = "distmat supplied")
    expect_identical(pc_dtwb@cluster, pc_dtwlb@cluster, info = "dtw_basic vs dtw_lb")
    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre@cluster, info = "no pam.pre but sparse")
    expect_identical(pc_dtwb@cluster, pc_dtwb_npampre2@cluster, info = "no pam.pre nor sparse")

    pc_dtwb <- reset_nondeterministic(pc_dtwb)
    pc_dtwb_npampre <- reset_nondeterministic(pc_dtwb_npampre)
    pc_dtwb_distmat <- reset_nondeterministic(pc_dtwb_distmat)
    pc_dtwlb <- reset_nondeterministic(pc_dtwlb)

    assign("pc_dtwb", pc_dtwb, persistent)
    assign("pc_dtwb_npampre", pc_dtwb_npampre, persistent)
    assign("pc_dtwb_distmat", pc_dtwb_distmat, persistent)
    assign("pc_dtwlb", pc_dtwlb, persistent)

    ## ---------------------------------------------------------- k-Shape
    pc_kshape <- tsclust(data_subset, type = "p", k = 4L,
                         distance = "sbd", centroid = "shape",
                         seed = 938)
    pc_kshape <- reset_nondeterministic(pc_kshape)
    assign("pc_kshape", pc_kshape, persistent)

    ## ---------------------------------------------------------- dba
    pc_dba <- tsclust(data_subset, type = "p", k = 4L,
                      distance = "dtw_basic", centroid = "dba",
                      args = tsclust_args(cent = list(max.iter = 15L)),
                      seed = 938)
    pc_dba <- reset_nondeterministic(pc_dba)
    assign("pc_dba", pc_dba, persistent)

    ## ---------------------------------------------------------- multivariate pam
    pc_mv_pam <- tsclust(data_multivariate, type = "p", k = 4L,
                         distance = "dtw_basic", centroid = "pam",
                         seed = 938)
    pc_mv_pam <- reset_nondeterministic(pc_mv_pam)
    assign("pc_mv_pam", pc_mv_pam, persistent)

    ## ---------------------------------------------------------- multivariate dba
    pc_mv_dba <- tsclust(data_multivariate, type = "p", k = 4L,
                         distance = "dtw_basic", centroid = "dba",
                         args = tsclust_args(cent = list(max.iter = 15L)),
                         seed = 938)
    pc_mv_dba <- reset_nondeterministic(pc_mv_dba)
    assign("pc_mv_dba", pc_mv_dba, persistent)

    ## ---------------------------------------------------------- sdtw
    pc_sdtw <- tsclust(data_multivariate, type = "p", k = 4L,
                       distance = "sdtw", centroid = "sdtw_cent",
                       seed = 938)
    pc_sdtw <- reset_nondeterministic(pc_sdtw)
    assign("pc_sdtw", pc_sdtw, persistent)
})

# ==================================================================================================
# TADPole
# ==================================================================================================

test_that("TADPole works as expected", {
    ## ---------------------------------------------------------- TADPole
    pc_tadp <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                       control = tadpole_control(dc = 1.5, window.size = 20L),
                       seed = 938)
    expect_s4_class(pc_tadp, "PartitionalTSClusters")
    pc_tadp <- reset_nondeterministic(pc_tadp)
    assign("pc_tadp", pc_tadp, persistent)

    ## ---------------------------------------------------------- TADPole with LBI
    pc_tadp_lbi <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                           control = tadpole_control(dc = 1.5, window.size = 20L, lb = "lbi"),
                           seed = 938)
    expect_s4_class(pc_tadp_lbi, "PartitionalTSClusters")
    pc_tadp_lbi <- reset_nondeterministic(pc_tadp_lbi)
    assign("pc_tadp_lbi", pc_tadp_lbi, persistent)

    ## ---------------------------------------------------------- TADPole with custom centroid
    pc_tadp_cent <- tsclust(data_reinterpolated_subset, type = "t", k = 4L,
                            preproc = zscore, centroid = shape_extraction,
                            control = tadpole_control(dc = 1.5, window.size = 20L),
                            seed = 938)
    expect_s4_class(pc_tadp_cent, "PartitionalTSClusters")
    pc_tadp_cent <- reset_nondeterministic(pc_tadp_cent)
    expect_identical(pc_tadp_cent@centroid, "shape_extraction")
    assign("pc_tadp_cent", pc_tadp_cent, persistent)

    pc_tadp_sdtwc <- tsclust(data_reinterpolated_subset, type = "t", k = 2L,
                             centroid = sdtw_cent,
                             control = tadpole_control(dc = 1.5, window.size = 20L),
                             seed = 938)
    expect_s4_class(pc_tadp_sdtwc, "PartitionalTSClusters")
    expect_identical(pc_tadp_sdtwc@centroid, "sdtw_cent")
})

# ==================================================================================================
# cluster reinitialization
# ==================================================================================================

## this case causes clusters to become empty
test_that("Cluster reinitialization in partitional tsclust works.", {
    expect_warning(
        pc_cr <- tsclust(data_reinterpolated, k = 20,
                         distance = "lbk", centroid = "mean",
                         seed = 31231,
                         control = partitional_control(iter.max = 10L, version = 1L),
                         args = tsclust_args(dist = list(window.size = 19L))),
        regexp = "converge"
    )
    expect_false(pc_cr@converged)
    pc_cr <- reset_nondeterministic(pc_cr)
    assign("pc_cr", pc_cr, persistent)

    ## test PAM too
    expect_warning(
        pc_cr <- tsclust(data_reinterpolated, k = 20,
                         distance = "lbk", centroid = "pam",
                         seed = 31231,
                         control = partitional_control(iter.max = 10L, version = 1L),
                         args = tsclust_args(dist = list(window.size = 19L))),
        regexp = "converge"
    )
    expect_false(pc_cr@converged)
})

# ==================================================================================================
# custom centroid
# ==================================================================================================

test_that("Operations with custom centroid complete successfully.", {
    ## ---------------------------------------------------------- with dots
    mycent <- function(x, cl_id, k, cent, cl_old, ...) {
        x_split <- split(x, cl_id)
        x_split <- lapply(x_split, function(xx) do.call(rbind, xx))
        new_cent <- lapply(x_split, colMeans)
        new_cent
    }

    cent_colMeans <- tsclust(data_matrix, k = 20L,
                             distance = "sbd", centroid = mycent, seed = 123)
    expect_s4_class(cent_colMeans, "PartitionalTSClusters")
    cent_colMeans <- reset_nondeterministic(cent_colMeans)

    ## ---------------------------------------------------------- without dots
    mycent <- function(x, cl_id, k, cent, cl_old) {
        x_split <- split(x, cl_id)
        x_split <- lapply(x_split, function(xx) do.call(rbind, xx))
        new_cent <- lapply(x_split, colMeans)
        new_cent
    }

    cent_colMeans_nd <- tsclust(data_matrix, k = 20L,
                                distance = "sbd", centroid = mycent, seed = 123)
    expect_s4_class(cent_colMeans_nd, "PartitionalTSClusters")
    cent_colMeans_nd <- reset_nondeterministic(cent_colMeans_nd)

    ## ---------------------------------------------------------- refs
    assign("cent_colMeans", cent_colMeans, persistent)
    assign("cent_colMeans_nd", cent_colMeans_nd, persistent)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
