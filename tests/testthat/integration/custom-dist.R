context("    Custom proxy distances and tsclust")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# test version without ellipsis
ndtw <- function(x, y, step.pattern = symmetric2,
                 window.type = "none", window.size = NULL,
                 open.end = FALSE, open.begin = FALSE)
{
    dtw::dtw(x, y, distance.only = TRUE,
             step.pattern = step.pattern,
             window.type = window.type,
             window.size = window.size,
             open.end = open.end,
             open.begin = open.begin)$normalizedDistance
}

if (!pr_DB$entry_exists("nDTW"))
    proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                           loop = TRUE, type = "metric", distance = TRUE,
                           description = "Normalized DTW with L1 norm")

# ==================================================================================================
# Registered with proxy
# ==================================================================================================

test_that("Calling tsclust after registering a custom distance works as expected.", {
    ## ---------------------------------------------------------- non-symmetric
    pc_ndtw <- tsclust(data_subset, k = 4L, distance = "nDTW", seed = 8319L,
                       control = partitional_control(version = 1L))
    pc_ndtw <- reset_nondeterministic(pc_ndtw)

    ## ---------------------------------------------------------- symmetric
    pc_ndtw_sym <- tsclust(data_subset, k = 4L, distance = "nDTW", seed = 8319L,
                           control = partitional_control(symmetric = TRUE,
                                                         version = 1L))
    pc_ndtw_sym <- reset_nondeterministic(pc_ndtw_sym)

    ## just for expect below
    pc_ndtw@control$symmetric <- TRUE
    pc_ndtw@call <- pc_ndtw_sym@call <- as.call(list("foo", bar = 1))

    expect_identical(pc_ndtw, pc_ndtw_sym)

    ## ---------------------------------------------------------- custom params
    pc_ndtw_par <- tsclust(data_subset, k = 4, distance = "nDTW", seed = 8319L,
                           args = tsclust_args(dist = list(window.type = "slantedband",
                                                           window.size = 18L,
                                                           open.begin = TRUE,
                                                           open.end = TRUE,
                                                           step.pattern = asymmetric)))

    pc_ndtw_par <- reset_nondeterministic(pc_ndtw_par)

    expect_equivalent(pc_ndtw_par@distmat,
                      proxy::dist(data_subset, data_subset, method = "nDTW",
                                  window.type = "slantedband", window.size = 18L,
                                  open.begin = TRUE, open.end = TRUE,
                                  step.pattern = asymmetric))

    assign("pc_ndtw", pc_ndtw, persistent)
    assign("pc_ndtw_sym", pc_ndtw_sym, persistent)
    assign("pc_ndtw_par", pc_ndtw_par, persistent)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
