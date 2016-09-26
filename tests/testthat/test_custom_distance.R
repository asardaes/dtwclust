context("Test custom proxy distances in dtwclust")

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

     pc_ndtw <- reset_nondeterministic(pc_ndtw)

     ## ---------------------------------------------------------- symmetric
     pc_ndtw_sym <- dtwclust(data_subset, k = 4, distance = "nDTW",
                             seed = 8319, control = list(symmetric = TRUE))

     pc_ndtw_sym <- reset_nondeterministic(pc_ndtw_sym)

     ## just for expect below
     pc_ndtw@control@symmetric <- TRUE
     pc_ndtw@call <- pc_ndtw_sym@call <- as.call(list("zas", a = 1))
     expect_identical(pc_ndtw, pc_ndtw_sym)

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

     skip_on_cran()

     expect_equal_to_reference(pc_ndtw, file_name(pc_ndtw), info = "nDTW")
     expect_equal_to_reference(pc_ndtw_sym, file_name(pc_ndtw_sym), info = "Symmetric nDTW")
     expect_equal_to_reference(pc_ndtw_par, file_name(pc_ndtw_par), info = "Params with nDTW")
})
