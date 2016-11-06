context("Test parallel")

# =================================================================================================
# Compare to sequential
# - This breaks 'CMD check' within RStudio, it just crashes when it creates the parallel workers
# - It still works if I manually build and check the project in the command line
# =================================================================================================

if (Sys.info()["user"] == "oso") {
    ## use all cores in local, change user if needed...
    num_workers <- parallel::detectCores()
} else {
    ## use only 2 in Travis CI
    num_workers <- 2L
}

## see https://github.com/hadley/testthat/issues/129
Sys.setenv("R_TESTS" = "")

test_that("Parallel computation gives the same results as sequential", {
    skip_on_cran()

    if (getOption("skip_par_tests", FALSE))
        skip("Parallel tests disabled explicitly.")

    cat("\n")

    require(doParallel)

    cl <- makeCluster(num_workers)
    invisible(clusterEvalQ(cl, library(dtwclust)))
    registerDoParallel(cl)

    ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
    test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

    stopCluster(cl)
    stopImplicitCluster()
    registerDoSEQ()

    rm(cl)
})

test_that("Parallel FORK computation gives the same results as sequential", {
    skip_on_cran()

    if (getOption("skip_par_tests", FALSE))
        skip("Parallel tests disabled explicitly.")

    skip_on_os(c("mac", "windows")) # mac in Travis CI is very slow

    ## Also test FORK in Linux
    cat("Test FORKs:\n")

    cl <- makeCluster(num_workers, "FORK")
    registerDoParallel(cl)

    ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
    test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

    stopCluster(cl)
    stopImplicitCluster()
    registerDoSEQ()
})
