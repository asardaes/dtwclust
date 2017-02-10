context("\tCustom proxy distances and dtwclust")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# Registered with proxy
# =================================================================================================

test_that("Calling dtwclust after registering a custom distance works as expected.", {
    ## ---------------------------------------------------------- nDTW
    ndtw <- function(x, y, ...) {
        dtw::dtw(x, y, distance.only = TRUE, ...)$normalizedDistance
    }

    if (!pr_DB$entry_exists("nDTW"))
        proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                               loop = TRUE, type = "metric", distance = TRUE,
                               description = "Normalized DTW with L1 norm")

    pc_ndtw <- dtwclust(data_subset, k = 4, distance = "nDTW", seed = 8319)
    pc_ndtw2 <- tsclust(data_subset, k = 4, distance = "nDTW", seed = 8319)

    pc_ndtw <- reset_nondeterministic(pc_ndtw)
    pc_ndtw2 <- reset_nondeterministic(pc_ndtw2)

    ## ---------------------------------------------------------- symmetric
    pc_ndtw_sym <- dtwclust(data_subset, k = 4, distance = "nDTW",
                            seed = 8319, control = list(symmetric = TRUE))

    pc_ndtw_sym2 <- tsclust(data_subset, k = 4, distance = "nDTW",
                            seed = 8319, control = partitional_control(symmetric = TRUE))

    pc_ndtw_sym <- reset_nondeterministic(pc_ndtw_sym)
    pc_ndtw_sym2 <- reset_nondeterministic(pc_ndtw_sym2)

    ## just for expect below
    pc_ndtw@control@symmetric <- TRUE
    pc_ndtw2@control$symmetric <- TRUE
    pc_ndtw@call <- pc_ndtw_sym@call <- as.call(list("zas", a = 1))
    pc_ndtw2@call <- pc_ndtw_sym2@call <- as.call(list("zas", a = 1))
    expect_identical(pc_ndtw, pc_ndtw_sym)
    expect_identical(pc_ndtw2, pc_ndtw_sym2)

    ## ---------------------------------------------------------- custom params
    pc_ndtw_par <- dtwclust(data_subset, k = 4, distance = "nDTW",
                            seed = 8319, control = list(window.size = 18L),
                            open.begin = TRUE, open.end = TRUE,
                            step.pattern = asymmetric)

    pc_ndtw_par <- reset_nondeterministic(pc_ndtw_par)

    expect_equal(pc_ndtw_par@distmat,
                 proxy::dist(data_subset, data_subset, method = "nDTW",
                             window.type = "slantedband", window.size = 18L,
                             open.begin = TRUE, open.end = TRUE,
                             step.pattern = asymmetric),
                 tolerance = 0,
                 check.attributes = FALSE)

    pc_ndtw_par2 <- tsclust(data_subset, k = 4, distance = "nDTW", seed = 8319,
                            args = tsclust_args(dist = list(window.size = 18L,
                                                            open.begin = TRUE,
                                                            open.end = TRUE,
                                                            step.pattern = asymmetric)))

    expect_identical(pc_ndtw_par@distmat, pc_ndtw_par2@distmat)

    assign("pc_ndtw", pc_ndtw, persistent)
    assign("pc_ndtw_sym", pc_ndtw_sym, persistent)
    assign("pc_ndtw_par", pc_ndtw_par, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
