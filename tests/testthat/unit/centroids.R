context("    Centroids")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

pt_ctrl <- partitional_control()
x <- data_reinterpolated_subset
k <- 2L
cl_id <- rep(c(1L, 2L), each = length(x) / 2L)
x_mv <- reinterpolate(data_multivariate, 205L)
expect_trace <- if (foreach::getDoParWorkers() > 1L) testthat::expect_silent else testthat::expect_output

# ==================================================================================================
# mean
# ==================================================================================================

test_that("Operations with mean centroid complete successfully.", {
    family <- new("tsclustFamily", control = pt_ctrl, allcent = "mean")

    ## ---------------------------------------------------------- univariate
    cent_mean <- family@allcent(x,
                                cl_id = cl_id,
                                k = k,
                                cent = x[c(1L,20L)],
                                cl_old = 0L)

    expect_identical(length(cent_mean), k)

    ## ---------------------------------------------------------- multivariate
    cent_mv_mean <- family@allcent(x_mv,
                                   cl_id = cl_id,
                                   k = k,
                                   cent = x_mv[c(1L,20L)],
                                   cl_old = 0L)

    expect_identical(length(cent_mv_mean), k)
    expect_identical(dim(cent_mv_mean[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_mean", cent_mean, persistent)
    assign("cent_mv_mean", cent_mv_mean, persistent)
})

# ==================================================================================================
# median
# ==================================================================================================

test_that("Operations with median centroid complete successfully.", {
    family <- new("tsclustFamily", control = pt_ctrl, allcent = "median")

    ## ---------------------------------------------------------- univariate
    cent_median <- family@allcent(x,
                                  cl_id = cl_id,
                                  k = k,
                                  cent = x[c(1L,20L)],
                                  cl_old = 0L)

    expect_identical(length(cent_median), k)

    ## ---------------------------------------------------------- multivariate
    cent_mv_median <- family@allcent(x_mv,
                                     cl_id = cl_id,
                                     k = k,
                                     cent = x_mv[c(1L,20L)],
                                     cl_old = 0L)

    expect_identical(length(cent_mv_median), k)
    expect_identical(dim(cent_mv_median[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_median", cent_median, persistent)
    assign("cent_mv_median", cent_mv_median, persistent)
})

# ==================================================================================================
# shape
# ==================================================================================================

test_that("Operations with shape centroid complete successfully.", {
    expect_error(shape_extraction(x_mv, x[[1L]]), "Dimension inconsistency")

    family <- new("tsclustFamily", control = pt_ctrl, allcent = "shape")

    ## ---------------------------------------------------------- univariate
    cent_shape <- family@allcent(x,
                                 cl_id = cl_id,
                                 k = k,
                                 cent = x[c(1L,20L)],
                                 cl_old = 0L)

    expect_identical(length(cent_shape), k)

    ## ---------------------------------------------------------- multivariate
    cent_mv_shape <- family@allcent(x_mv,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x_mv[c(1L,20L)],
                                    cl_old = 0L)

    expect_identical(length(cent_mv_shape), k)
    expect_identical(dim(cent_mv_shape[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_shape", cent_shape, persistent)
    assign("cent_mv_shape", cent_mv_shape, persistent)

    ## ---------------------------------------------------------- directly
    zero_shape_centroid <- shape_extraction(matrix(0, 2L, 10L))
    expect_identical(sum(zero_shape_centroid), 0)
    zero_shape_centroid_reference <- shape_extraction(matrix(0, 2L, 10L), x[[1L]])
    expect_identical(zero_shape_centroid_reference, x[[1L]])

    expect_error(shape_extraction(data_subset, data_multivariate[[1L]]), regexp = "dimension")
})

# ==================================================================================================
# pam
# ==================================================================================================

test_that("Operations with pam centroid complete successfully.", {
    ## ---------------------------------------------------------- univariate without distmat
    dm <- dtwclust:::Distmat$new(series = x,
                                 control = pt_ctrl,
                                 distance = "sbd")
    pt_ctrl$distmat <- dm
    family <- new("tsclustFamily",
                  control = pt_ctrl,
                  dist = "sbd",
                  allcent = "pam")

    cent_pam <- family@allcent(x,
                               cl_id = cl_id,
                               k = k,
                               cent = x[c(1L,20L)],
                               cl_old = 0L)

    expect_identical(length(cent_pam), k)
    expect_null(as.list(environment(family@allcent))$distmat)

    ## ---------------------------------------------------------- univariate with R distmat
    pt_ctrl$distmat <- proxy::dist(x, method = "sbd")
    family <- new("tsclustFamily",
                  control = pt_ctrl,
                  dist = "sbd",
                  allcent = "pam")

    cent_pam_distmat <- family@allcent(x,
                                       cl_id = cl_id,
                                       k = k,
                                       cent = x[c(1L,20L)],
                                       cl_old = 0L)

    expect_identical(length(cent_pam_distmat), k)
    expect_false(is.null(environment(family@allcent)$control$distmat))
    expect_identical(cent_pam, cent_pam_distmat)

    ## ---------------------------------------------------------- multivariate with R distmat
    pt_ctrl$distmat <- proxy::dist(x_mv, method = "dtw_basic", window.size = 18L)
    family <- new("tsclustFamily",
                  control = pt_ctrl,
                  dist = "dtw_basic",
                  allcent = "pam")

    cent_mv_pam <- family@allcent(x_mv,
                                  cl_id = cl_id,
                                  k = k,
                                  cent = x_mv[c(1L,20L)],
                                  cl_old = 0L)

    expect_identical(length(cent_mv_pam), k)
    expect_identical(dim(cent_mv_pam[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- sparse symmetric
    pt_ctrl$symmetric <- TRUE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = pt_ctrl,
                                       distance = "sbd",
                                       dist_args = list())
    pt_ctrl$distmat <- dm
    family <- new("tsclustFamily",
                  control = pt_ctrl,
                  dist = "sbd",
                  allcent = "pam")

    expect_equal(cent_pam,
                 family@allcent(x,
                                cl_id = cl_id,
                                k = k,
                                cent = x[c(1L,20L)],
                                cl_old = 0L),
                 check.attributes = FALSE)

    ## ---------------------------------------------------------- sparse non-symmetric
    pt_ctrl$symmetric <- FALSE
    dm <- dtwclust:::SparseDistmat$new(series = x_mv,
                                       control = pt_ctrl,
                                       distance = "dtw_basic",
                                       dist_args = list(window.size = 18L))
    pt_ctrl$distmat <- dm
    family <- new("tsclustFamily",
                  control = pt_ctrl,
                  dist = "dtw_basic",
                  allcent = "pam")

    expect_equal(cent_mv_pam,
                 family@allcent(x_mv,
                                cl_id = cl_id,
                                k = k,
                                cent = x_mv[c(1L,20L)],
                                cl_old = 0L),
                 check.attributes = FALSE)

    ## ---------------------------------------------------------- refs
    assign("cent_pam", cent_pam, persistent)
    assign("cent_mv_pam", cent_mv_pam, persistent)

    ## ---------------------------------------------------------- standalone
    expect_error(pam_cent(x))

    dm <- proxy::dist(x[6L:10L], method = "dtw_basic", window.size = 15L)
    pam_cent_no_distmat <- pam_cent(x, "dtw_basic", 6L:10L, window.size = 15L)
    expect_identical(attr(pam_cent_no_distmat, "series_id"), 7L)
    pam_cent_with_distmat <- pam_cent(x[6L:10L], distmat = dm)
    expect_identical(attr(pam_cent_with_distmat, "series_id"), 2L)
    expect_equal(pam_cent_with_distmat, pam_cent_no_distmat, check.attributes = FALSE)
})

# ==================================================================================================
# dba
# ==================================================================================================

test_that("Operations with dba centroid complete successfully.", {
    expect_error(dba(x[1L:5L], max.iter = 0L), "positive")
    expect_error(dba(x[1L:5L], step.pattern = dtw::asymmetric), "step.pattern")

    input_order1 <- dba(x[1L:5L], x[[1L]], step.pattern = dtw::symmetric1)
    input_order2 <- dba(x[5L:1L], x[[1L]], step.pattern = dtw::symmetric1)
    expect_equal(input_order1, input_order2)

    one_delta <-  dba(x[1L:5L], x[[1L]], delta = 0.01)
    many_delta <-  dba(x[1L:5L], x[[1L]], delta = c(0.01, 10))
    expect_equal(one_delta, many_delta)

    family <- new("tsclustFamily", control = pt_ctrl, allcent = "dba")

    ## ---------------------------------------------------------- univariate

    expect_output(
        dba(x[1L:5L], centroid = x[[1L]], max.iter = 20L, trace = TRUE),
        "Converged!"
    )

    expect_trace(
        cent_dba <- family@allcent(x,
                                   cl_id = cl_id,
                                   k = k,
                                   cent = x[c(1L,20L)],
                                   cl_old = 0L,
                                   window.size = 18L,
                                   max.iter = 15L,
                                   trace = TRUE)
    )

    expect_identical(length(cent_dba), k)

    ## ---------------------------------------------------------- multivariate
    expect_trace(
        cent_mv_dba <- family@allcent(x_mv,
                                      cl_id = cl_id,
                                      k = k,
                                      cent = x_mv[c(1L,20L)],
                                      cl_old = 0L,
                                      window.size = 18L,
                                      max.iter = 15L,
                                      norm = "L2",
                                      trace = TRUE)
    )

    expect_identical(length(cent_mv_dba), k)
    expect_identical(dim(cent_mv_dba[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- multivariate v2
    expect_trace(
        cent_mv_dba_bys <- family@allcent(x_mv,
                                          cl_id = cl_id,
                                          k = k,
                                          cent = x_mv[c(1L,20L)],
                                          cl_old = 0L,
                                          mv.ver = "by-s",
                                          trace = TRUE)
    )

    expect_identical(length(cent_mv_dba_bys), k)
    expect_identical(dim(cent_mv_dba_bys[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_dba", cent_dba, persistent)
    assign("cent_mv_dba", cent_mv_dba, persistent)
    assign("cent_mv_dba_bys", cent_mv_dba_bys, persistent)
})

# ==================================================================================================
# sdtw_cent
# ==================================================================================================

test_that("Operations with sdtw_cent centroid complete successfully.", {
    expect_error(sdtw_cent(x, gamma = -1), "gamma.*positive")
    expect_error(sdtw_cent(x, weights = 1), "weights.*length")

    family <- new("tsclustFamily", control = pt_ctrl, allcent = "sdtw_cent")

    ## ---------------------------------------------------------- univariate

    cent_sdtwc <- family@allcent(x,
                                 cl_id = cl_id,
                                 k = k,
                                 cent = x[c(1L,20L)],
                                 cl_old = 0L)

    expect_identical(length(cent_sdtwc), k)

    ## ---------------------------------------------------------- multivariate
    cent_mv_sdtwc <- family@allcent(x_mv,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x_mv[c(1L,20L)],
                                    cl_old = 0L)

    expect_identical(length(cent_mv_sdtwc), k)
    expect_identical(dim(cent_mv_sdtwc[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_sdtwc", cent_sdtwc, persistent)
    assign("cent_mv_sdtwc", cent_mv_sdtwc, persistent)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
