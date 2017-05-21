context("\tCentroids")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

ctrl <- new("dtwclustControl", window.size = 18L)
pt_ctrl <- partitional_control()
x <- data_reinterpolated_subset
k <- 2L
cl_id <- rep(c(1L, 2L), each = length(x) / 2L)
x_mv <- reinterpolate(data_multivariate, 205L)

# =================================================================================================
# mean
# =================================================================================================

test_that("Operations with mean centroid complete successfully.", {
    ## ---------------------------------------------------------- univariate
    family <- new("dtwclustFamily",
                  control = ctrl,
                  allcent = "mean")

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

    ## ---------------------------------------------------------- tsclustFamily
    family <- new("tsclustFamily", control = pt_ctrl, allcent = "mean")

    expect_identical(cent_mean,
                     family@allcent(x,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x[c(1L,20L)],
                                    cl_old = 0L))

    expect_identical(cent_mv_mean,
                     family@allcent(x_mv,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x_mv[c(1L,20L)],
                                    cl_old = 0L))
})

# =================================================================================================
# median
# =================================================================================================

test_that("Operations with median centroid complete successfully.", {
    ## ---------------------------------------------------------- univariate
    family <- new("dtwclustFamily",
                  control = ctrl,
                  allcent = "median")

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

    ## ---------------------------------------------------------- tsclustFamily
    family <- new("tsclustFamily", control = pt_ctrl, allcent = "median")

    expect_identical(cent_median,
                     family@allcent(x,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x[c(1L,20L)],
                                    cl_old = 0L))

    expect_identical(cent_mv_median,
                     family@allcent(x_mv,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x_mv[c(1L,20L)],
                                    cl_old = 0L))
})

# =================================================================================================
# shape
# =================================================================================================

test_that("Operations with shape centroid complete successfully.", {
    ## ---------------------------------------------------------- univariate
    family <- new("dtwclustFamily",
                  control = ctrl,
                  allcent = "shape")

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

    ## ---------------------------------------------------------- tsclustFamily
    family <- new("tsclustFamily", control = pt_ctrl, allcent = "shape")

    expect_identical(cent_shape,
                     family@allcent(x,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x[c(1L,20L)],
                                    cl_old = 0L))

    expect_identical(cent_mv_shape,
                     family@allcent(x_mv,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x_mv[c(1L,20L)],
                                    cl_old = 0L))
})

# =================================================================================================
# pam
# =================================================================================================

test_that("Operations with pam centroid complete successfully.", {
    ## ---------------------------------------------------------- univariate without distmat
    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "sbd",
                  allcent = "pam")

    cent_pam <- family@allcent(x,
                               cl_id = cl_id,
                               k = k,
                               cent = x[c(1L,20L)],
                               cl_old = 0L)

    expect_identical(length(cent_pam), k)

    expect_null(as.list(environment(family@allcent))$distmat)

    ## ---------------------------------------------------------- univariate with distmat
    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "sbd",
                  allcent = "pam",
                  distmat = proxy::dist(x, method = "sbd"))

    cent_pam_distmat <- family@allcent(x,
                                       cl_id = cl_id,
                                       k = k,
                                       cent = x[c(1L,20L)],
                                       cl_old = 0L)

    expect_identical(length(cent_pam_distmat), k)

    expect_false(is.null(as.list(environment(family@allcent))$distmat))

    expect_identical(cent_pam, cent_pam_distmat)

    ## ---------------------------------------------------------- multivariate
    family <- new("dtwclustFamily",
                  control = ctrl,
                  dist = "dtw_basic",
                  allcent = "pam")

    cent_mv_pam <- family@allcent(x_mv,
                                  cl_id = cl_id,
                                  k = k,
                                  cent = x_mv[c(1L,20L)],
                                  cl_old = 0L)

    expect_identical(length(cent_mv_pam), k)

    expect_identical(dim(cent_mv_pam[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_pam", cent_pam, persistent)
    assign("cent_mv_pam", cent_mv_pam, persistent)

    ## ---------------------------------------------------------- tsclustFamily

    ## with distmat
    pt_ctrl$distmat <- proxy::dist(x, method = "sbd")
    family <- new("tsclustFamily",
                  control = pt_ctrl,
                  dist = "sbd",
                  allcent = "pam")
    expect_false(is.null(as.list(environment(family@allcent))$control$distmat))

    expect_equal(cent_pam_distmat,
                 family@allcent(x,
                                cl_id = cl_id,
                                k = k,
                                cent = x[c(1L,20L)],
                                cl_old = 0L),
                 check.attributes = FALSE)

    ## multivariate
    pt_ctrl$distmat <- proxy::dist(x_mv, method = "dtw_basic", window.size = 18L)
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

    ## sparse symmetric
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

    ## sparse non-symmetric
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
})

# =================================================================================================
# dba
# =================================================================================================

test_that("Operations with dba centroid complete successfully.", {
    ## ---------------------------------------------------------- univariate
    family <- new("dtwclustFamily",
                  control = ctrl,
                  allcent = "dba")

    cent_dba <- family@allcent(x,
                               cl_id = cl_id,
                               k = k,
                               cent = x[c(1L,20L)],
                               cl_old = 0L)

    expect_identical(length(cent_dba), k)

    ## ---------------------------------------------------------- multivariate
    ctrl@norm <- "L2"
    family2 <- new("dtwclustFamily",
                   control = ctrl,
                   allcent = "dba")

    cent_mv_dba <- family2@allcent(x_mv,
                                   cl_id = cl_id,
                                   k = k,
                                   cent = x_mv[c(1L,20L)],
                                   cl_old = 0L)

    expect_identical(length(cent_mv_dba), k)

    expect_identical(dim(cent_mv_dba[[1L]]), dim(x_mv[[1L]]))

    ## ---------------------------------------------------------- refs
    assign("cent_dba", cent_dba, persistent)
    assign("cent_mv_dba", cent_mv_dba, persistent)

    ## ---------------------------------------------------------- tsclustFamily
    family <- new("tsclustFamily", control = pt_ctrl, allcent = "dba")

    expect_identical(cent_dba,
                     family@allcent(x,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x[c(1L,20L)],
                                    cl_old = 0L,
                                    window.size = 18L,
                                    max.iter = ctrl@dba.iter))

    expect_identical(cent_mv_dba,
                     family@allcent(x_mv,
                                    cl_id = cl_id,
                                    k = k,
                                    cent = x_mv[c(1L,20L)],
                                    cl_old = 0L,
                                    window.size = 18L,
                                    max.iter = ctrl@dba.iter,
                                    norm = "L2"))
})

test_that("Second version of DBA works as expected.", {
    dba2_uv <- DBA(data_subset[1L:5L], centroid = data_subset[[1L]], mv.ver = "by-s")
    expect_true(is.null(dim(dba2_uv)))

    dba2_mv <- DBA(data_multivariate[1L:5L], centroid = data_multivariate[[1L]], mv.ver = "by-s")
    expect_identical(dim(dba2_mv), dim(data_multivariate[[1L]]))
    expect_identical(dimnames(dba2_mv), dimnames(data_multivariate[[1L]]))

    assign("dba2_uv", dba2_uv, persistent)
    assign("dba2_mv", dba2_mv, persistent)
})

# =================================================================================================
# custom
# =================================================================================================

test_that("Operations with custom centroid complete successfully.", {
    ## ---------------------------------------------------------- with dots
    mycent <- function(x, cl_id, k, cent, cl_old, ...) {
        x_split <- split(x, cl_id)

        x_split <- lapply(x_split, function(xx) do.call(rbind, xx))

        new_cent <- lapply(x_split, colMeans)

        new_cent
    }

    cent_colMeans <- dtwclust(data_matrix, type = "partitional", k = 20,
                              distance = "sbd", centroid = mycent,
                              preproc = NULL, control = ctrl, seed = 123)

    cent_colMeans <- reset_nondeterministic(cent_colMeans)

    ## ---------------------------------------------------------- without dots
    mycent <- function(x, cl_id, k, cent, cl_old) {
        x_split <- split(x, cl_id)

        x_split <- lapply(x_split, function(xx) do.call(rbind, xx))

        new_cent <- lapply(x_split, colMeans)

        new_cent
    }

    cent_colMeans_nd <- dtwclust(data_matrix, type = "partitional", k = 20L,
                                 distance = "sbd", centroid = mycent,
                                 preproc = NULL, control = ctrl, seed = 123)

    cent_colMeans_nd <- reset_nondeterministic(cent_colMeans)

    ## ---------------------------------------------------------- refs
    assign("cent_colMeans", cent_colMeans, persistent)
    assign("cent_colMeans_nd", cent_colMeans_nd, persistent)

    ## ---------------------------------------------------------- tsclustFamily
    pc <- tsclust(data_matrix, k = 20L,
                  distance = "sbd", centroid = mycent,
                  seed = 123,
                  args = tsclust_args(dist = list(window.size = 18L)))

    expect_identical(cent_colMeans_nd@centroids, pc@centroids)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
