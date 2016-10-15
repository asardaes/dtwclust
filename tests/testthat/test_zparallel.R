context("Test parallel")

# =================================================================================================
# Compare to sequential
# - This breaks 'CMD check' within RStudio, it just crashes when it creates the parallel workers
# - It still works if I manually build and check the project in the command line
# =================================================================================================

test_that("Parallel computation gives the same results as sequential", {
    skip_on_cran()

    if (getOption("skip_par_tests", FALSE))
        skip("Parallel tests disabled explicitly.")

    ## see https://github.com/hadley/testthat/issues/129
    Sys.setenv("R_TESTS" = "")

    cat("\n")

    require(doParallel)

    if (identical(Sys.getenv("NOT_CRAN"), ""))
        num_workers <- detectCores()
    else
        num_workers <- 2L

    cl <- makeCluster(num_workers)
    invisible(clusterEvalQ(cl, library(dtwclust)))
    registerDoParallel(cl)

    ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
    test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

    stopCluster(cl)
    stopImplicitCluster()
    registerDoSEQ()

    skip_on_os("mac")

    ## Also test FORK in Unix
    if (tolower(Sys.info()[["sysname"]]) != "windows") {
        cat("Test FORKs:\n")

        rm(cl)
        cl <- makeCluster(num_workers, "FORK")
        registerDoParallel(cl)

        ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
        test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

        stopCluster(cl)
        stopImplicitCluster()
        registerDoSEQ()
    }
})
