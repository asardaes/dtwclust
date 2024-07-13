# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# preproc
# ==================================================================================================

test_that("Preprocessing functions are called by tsclust without errors.", {
    pc0 <- tsclust(data_subset, k = 4L,
                   distance = "sbd", centroid = "shape",
                   keep.attributes = TRUE,
                   seed = 1899)

    preproc_dots <- function(x, ...) {
        zscore(x, ...)
    }

    preproc_nodots <- function(x, keep.attributes = FALSE) {
        zscore(x, keep.attributes = keep.attributes)
    }

    pc1 <- tsclust(data_subset, k = 4L,
                   distance = "sbd", centroid = "shape",
                   preproc = preproc_dots, keep.attributes = TRUE,
                   seed = 1899)

    pc2 <- tsclust(data_subset, k = 4L,
                   distance = "sbd", centroid = "shape",
                   preproc = preproc_nodots, keep.attributes = TRUE,
                   seed = 1899)

    pc0 <- reset_nondeterministic(pc0)
    pc1 <- reset_nondeterministic(pc1)
    pc2 <- reset_nondeterministic(pc2)

    pc0@call <- pc1@call <- pc2@call <- call("foo", x = 1)
    pc1@preproc <- pc2@preproc <- pc0@preproc
    # this is not identical with MKL? why?
    expect_equal(pc0, pc1)
    expect_equal(pc0, pc2)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
