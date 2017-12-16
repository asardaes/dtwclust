context("    Families and proxy distances")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

ts_ctrl <- partitional_control()
x <- data_reinterpolated_subset
centroids <- x[c(1L, 15L)]
attr(centroids, "id_cent") <- c(1L, 15L)
window.size <- 18L

# ==================================================================================================
# lbk
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and lbk give expected results", {
    distmat <- proxy::dist(x, x, method = "lbk", window.size = window.size)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "lbk")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "lbk")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_lbk", whole_distmat, persistent)

    ## ---------------------------------------------------------- tsclustFamily, sparse distmat
    ts_ctrl$symmetric <- FALSE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = ts_ctrl,
                                       distance = "lbk",
                                       dist_args = list(window.size = window.size),
                                       id_cent = attr(centroids, "id_cent"))
    ts_ctrl$distmat <- dm

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "lbk")

    sub_distmat <- base::as.matrix(family@dist(x, centroids))
    whole_distmat <- base::as.matrix(family@dist(x))

    expect_equal(whole_distmat, base::as.matrix(distmat), info = "Whole, sparse distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, base::as.matrix(distmat[ , c(1L, 15L), drop = FALSE]),
                 info = "Sub, sparse distmat", tolerance = 0, check.attributes = FALSE)
})

# ==================================================================================================
# lbi
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and lbi give expected results", {
    distmat <- proxy::dist(x, x, method = "lbi", window.size = window.size)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "lbi")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "lbi")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_lbi", whole_distmat, persistent)

    ## ---------------------------------------------------------- tsclustFamily, sparse distmat
    ts_ctrl$symmetric <- FALSE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = ts_ctrl,
                                       distance = "lbi",
                                       dist_args = list(window.size = window.size),
                                       id_cent = attr(centroids, "id_cent"))
    ts_ctrl$distmat <- dm

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "lbi")

    sub_distmat <- base::as.matrix(family@dist(x, centroids))
    whole_distmat <- base::as.matrix(family@dist(x))

    expect_equal(whole_distmat, base::as.matrix(distmat), info = "Whole, sparse distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, base::as.matrix(distmat[ , c(1L, 15L), drop = FALSE]),
                 info = "Sub, sparse distmat", tolerance = 0, check.attributes = FALSE)
})

# ==================================================================================================
# sbd
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and sbd give expected results", {
    distmat <- proxy::dist(x, x, method = "sbd")

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "sbd")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "sbd")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_sbd", whole_distmat, persistent)

    ## ---------------------------------------------------------- tsclustFamily, sparse distmat
    ts_ctrl$symmetric <- TRUE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = ts_ctrl,
                                       distance = "sbd",
                                       dist_args = list(),
                                       id_cent = attr(centroids, "id_cent"))
    ts_ctrl$distmat <- dm

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "sbd")

    sub_distmat <- base::as.matrix(family@dist(x, centroids))
    whole_distmat <- base::as.matrix(family@dist(x))

    expect_equal(whole_distmat, base::as.matrix(distmat), info = "Whole, sparse distmat",
                 check.attributes = FALSE)

    expect_equal(sub_distmat, base::as.matrix(distmat[ , c(1L, 15L), drop = FALSE]),
                 info = "Sub, sparse distmat", check.attributes = FALSE)
})

# ==================================================================================================
# dtw_lb
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and dtw_lb give expected results", {
    distmat <- proxy::dist(x, x, method = "dtw_lb", window.size = window.size)
    sdm <- proxy::dist(x, centroids, method = "dtw_lb", window.size = window.size)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw_lb")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE)
    class(pdist) <- NULL

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, sdm, tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw_lb")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    sdm <- distmat[ , c(1L, 15L), drop = FALSE]
    dimnames(sdm) <- NULL
    expect_equal(sub_distmat, sdm, info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_dtwlb", whole_distmat, persistent)
})

# ==================================================================================================
# dtw
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and dtw give expected results", {
    distmat <- proxy::dist(x, x, method = "dtw", window.type = "slantedband", window.size = window.size)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_dtw", whole_distmat, persistent)
})

# ==================================================================================================
# dtw2
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and dtw2 give expected results", {
    distmat <- proxy::dist(x, x, method = "dtw2", window.type = "slantedband", window.size = window.size)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw2")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw2")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_dtw2", whole_distmat, persistent)
})

# ==================================================================================================
# dtw_basic
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and dtw_basic give expected results", {
    distmat <- proxy::dist(x, x, method = "dtw_basic", window.size = window.size)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw_basic")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw_basic")

    whole_distmat <- family@dist(x, window.size = window.size)
    sub_distmat <- family@dist(x, centroids, window.size = window.size)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_dtwb", whole_distmat, persistent)

    ## ---------------------------------------------------------- tsclustFamily, sparse distmat
    ts_ctrl$symmetric <- FALSE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = ts_ctrl,
                                       distance = "dtw_basic",
                                       dist_args = list(window.size = window.size),
                                       id_cent = attr(centroids, "id_cent"))
    ts_ctrl$distmat <- dm

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "dtw_basic")

    sub_distmat <- base::as.matrix(family@dist(x, centroids))
    whole_distmat <- base::as.matrix(family@dist(x))

    expect_equal(whole_distmat, base::as.matrix(distmat), info = "Whole, sparse distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, base::as.matrix(distmat[ , c(1L, 15L), drop = FALSE]),
                 info = "Sub, sparse distmat", tolerance = 0, check.attributes = FALSE)
})

# ==================================================================================================
# gak
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and gak give expected results", {
    distmat <- proxy::dist(x, x, method = "gak",
                           window.size = window.size,
                           sigma = 100)

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "gak")

    whole_distmat <- family@dist(x, window.size = window.size, sigma = 100)
    sub_distmat <- family@dist(x, centroids, window.size = window.size, sigma = 100)
    pdist <- family@dist(x, window.size = window.size, pairwise = TRUE, sigma = 100)
    class(pdist) <- NULL
    class(sub_distmat) <- "matrix"

    expect_equal(pdist, rep(0, length(pdist)), info = "Pairwise",
                 check.attributes = FALSE)

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "gak")

    whole_distmat <- family@dist(x, window.size = window.size, sigma = 100)
    sub_distmat <- family@dist(x, centroids, window.size = window.size, sigma = 100)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_gak", whole_distmat, persistent)

    ## ---------------------------------------------------------- tsclustFamily, sparse distmat
    ts_ctrl$symmetric <- FALSE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = ts_ctrl,
                                       distance = "gak",
                                       dist_args = list(window.size = window.size,
                                                        sigma = 100),
                                       id_cent = attr(centroids, "id_cent"))
    ts_ctrl$distmat <- dm

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "gak")

    sub_distmat <- base::as.matrix(family@dist(x, centroids))
    whole_distmat <- base::as.matrix(family@dist(x))

    expect_equal(whole_distmat, base::as.matrix(distmat), info = "Whole, sparse distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, base::as.matrix(distmat[ , c(1L, 15L), drop = FALSE]),
                 info = "Sub, sparse distmat", tolerance = 0, check.attributes = FALSE)
})

# ==================================================================================================
# sdtw
# ==================================================================================================

test_that("Operations with tsclustFamily@dist and sdtw give expected results", {
    distmat <- proxy::dist(x, x, method = "sdtw")

    ## ---------------------------------------------------------- tsclustFamily, no distmat
    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "sdtw")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)
    pdist <- family@dist(x, pairwise = TRUE)
    class(sub_distmat) <- "matrix"

    # different because sdtw(x,x) != 0
    expect_identical(length(pdist), length(x), info = "Pairwise")

    expect_equal(whole_distmat, distmat, info = "Whole, NULL distmat",
                 check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, NULL distmat",
                 check.attributes = FALSE)

    ## ---------------------------------------------------------- tsclustFamily, with distmat
    ts_ctrl$distmat <- dtwclust:::Distmat$new(distmat = distmat,
                                              id_cent = attr(centroids, "id_cent"))

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "sdtw")

    whole_distmat <- family@dist(x)
    sub_distmat <- family@dist(x, centroids)

    expect_equal(whole_distmat, distmat, info = "Whole, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    expect_equal(sub_distmat, distmat[ , c(1L, 15L), drop = FALSE], info = "Sub, with distmat",
                 tolerance = 0, check.attributes = FALSE)

    ## ---------------------------------------------------------- ref
    assign("distmat_sdtw", whole_distmat, persistent)

    ## ---------------------------------------------------------- tsclustFamily, sparse distmat
    ts_ctrl$symmetric <- FALSE
    dm <- dtwclust:::SparseDistmat$new(series = x,
                                       control = ts_ctrl,
                                       distance = "sdtw",
                                       id_cent = attr(centroids, "id_cent"))
    ts_ctrl$distmat <- dm

    family <- new("tsclustFamily",
                  control = ts_ctrl,
                  dist = "sdtw")

    sub_distmat <- base::as.matrix(family@dist(x, centroids))
    whole_distmat <- base::as.matrix(family@dist(x))

    expect_equal(whole_distmat, base::as.matrix(distmat), info = "Whole, sparse distmat",
                 check.attributes = FALSE)

    expect_equal(sub_distmat, base::as.matrix(distmat[ , c(1L, 15L), drop = FALSE]),
                 info = "Sub, sparse distmat", check.attributes = FALSE)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
