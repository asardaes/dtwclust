context("Test partitional algorithms")

# =================================================================================================
# invalid combinations
# =================================================================================================

test_that("Invalid combinations in partitional clustering are detected.", {
     expect_error(dtwclust(data, k = 101L), "more clusters")

     expect_error(dtwclust(data, distance = "lbk"), "different length")
     expect_error(dtwclust(data, distance = "lbi"), "different length")
     expect_error(dtwclust(data, distance = "dtw_lb"), "different length")

     expect_error(dtwclust(data, centroid = "mean"), "different length")
     expect_error(dtwclust(data, centroid = "median"), "different length")
     expect_error(dtwclust(data, centroid = "fcm"), "arg")

     expect_error(dtwclust(data, preproc = "zscore"), "preprocessing")
})

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

     skip_on_cran()

     expect_equal_to_reference(pc_k, file_name(pc_k))
     expect_equal_to_reference(pc_rep, file_name(pc_rep))
     expect_equal_to_reference(pc_krep, file_name(pc_krep))
})

# =================================================================================================
# partitional algorithms
# =================================================================================================

test_that("Partitional clustering works as expected.", {
     ## ---------------------------------------------------------- dtw
     pc_dtwb <- dtwclust(data_reinterpolated[1L:20L], type = "p", k = 4L,
                         distance = "dtw_basic", centroid = "pam",
                         seed = 938, control = list(window.size = 20L))

     pc_dtwb_npampre <- dtwclust(data_reinterpolated[1L:20L], type = "p", k = 4L,
                                 distance = "dtw_basic", centroid = "pam",
                                 seed = 938, control = list(window.size = 20L,
                                                            pam.precompute = FALSE))

     pc_dtwb_distmat <- dtwclust(data_reinterpolated[1L:20L], type = "p", k = 4L,
                                 distance = "dtw_basic", centroid = "pam",
                                 seed = 938, control = list(window.size = 20L),
                                 distmat = pc_dtwb@distmat)

     pc_dtwlb <- dtwclust(data_reinterpolated[1L:20L], type = "p", k = 4L,
                          distance = "dtw_lb", centroid = "pam",
                          seed = 938, control = list(window.size = 20L,
                                                     pam.precompute = FALSE))

     expect_identical(pc_dtwb@cluster, pc_dtwb_npampre@cluster, info = "pam.pre vs no pam.pre")
     expect_identical(pc_dtwb@cluster, pc_dtwb_distmat@cluster, info = "distmat supplied")
     expect_identical(pc_dtwb@cluster, pc_dtwlb@cluster, info = "dtw_basic vs dtw_lb")

     skip_on_cran()

     pc_dtwb <- reset_nondeterministic(pc_dtwb)
     pc_dtwb_npampre <- reset_nondeterministic(pc_dtwb_npampre)
     pc_dtwb_distmat <- reset_nondeterministic(pc_dtwb_distmat)
     pc_dtwlb <- reset_nondeterministic(pc_dtwlb)

     expect_equal_to_reference(pc_dtwb, file_name(pc_dtwb))
     expect_equal_to_reference(pc_dtwb_npampre, file_name(pc_dtwb_npampre))
     expect_equal_to_reference(pc_dtwb_distmat, file_name(pc_dtwb_distmat))
     expect_equal_to_reference(pc_dtwlb, file_name(pc_dtwlb))

     ## ---------------------------------------------------------- k-Shape
     pc_kshape <- dtwclust(data_subset, type = "p", k = 4L,
                           distance = "sbd", centroid = "shape",
                           seed = 938)

     pc_kshape <- reset_nondeterministic(pc_kshape)

     expect_equal_to_reference(pc_kshape, file_name(pc_kshape))

     ## ---------------------------------------------------------- dba
     pc_dba <- dtwclust(data_subset, type = "p", k = 4L,
                        distance = "dtw_basic", centroid = "dba",
                        seed = 938)

     pc_dba <- reset_nondeterministic(pc_dba)

     expect_equal_to_reference(pc_dba, file_name(pc_dba))

     ## ---------------------------------------------------------- multivariate pam
     pc_mv_pam <- dtwclust(data_multivariate[1L:20L], type = "p", k = 4L,
                           distance = "dtw_basic", centroid = "pam",
                           seed = 938)

     pc_mv_pam <- reset_nondeterministic(pc_mv_pam)

     expect_equal_to_reference(pc_mv_pam, file_name(pc_mv_pam))

     ## ---------------------------------------------------------- multivariate dba
     pc_mv_dba <- dtwclust(data_multivariate[1L:20L], type = "p", k = 4L,
                           distance = "dtw_basic", centroid = "dba",
                           seed = 938)

     pc_mv_dba <- reset_nondeterministic(pc_mv_dba)

     expect_equal_to_reference(pc_mv_dba, file_name(pc_mv_dba))
})

# =================================================================================================
# TADPole
# =================================================================================================

test_that("TADPole works as expected", {
     expect_error(dtwclust(data, type = "t", k = 20, control = list(window.size = 20L)), "dc")
     expect_error(dtwclust(data, type = "t", k = 20, dc = 1.5), "window.size")
     expect_error(dtwclust(data, type = "t", k = 20, dc = 1.5, control = list(window.size = 20L)),
                  "same length")

     skip_on_cran()

     ## ---------------------------------------------------------- TADPole
     pc_tadp <- dtwclust(data_reinterpolated[1L:20L], type = "t", k = 4L,
                         dc = 1.5, control = list(window.size = 20L),
                         seed = 938)

     pc_tadp <- reset_nondeterministic(pc_tadp)

     expect_equal_to_reference(pc_tadp, file_name(pc_tadp))

     ## ---------------------------------------------------------- TADPole with custom centroid
     pc_tadp_cent <- dtwclust(data_reinterpolated[1L:20L], type = "t", k = 4L,
                              preproc = zscore, centroid = shape_extraction,
                              dc = 1.5, control = list(window.size = 20L),
                              seed = 938)

     pc_tadp_cent <- reset_nondeterministic(pc_tadp_cent)

     expect_equal_to_reference(pc_tadp_cent, file_name(pc_tadp_cent))
})
