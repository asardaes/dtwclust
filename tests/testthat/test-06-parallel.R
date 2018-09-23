context("Parallel tests")

# =================================================================================================
# run all tests with a parallel backend
# =================================================================================================

chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CHECK/Travis/AppVeyor
    num_workers <- 2L
} else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
}

test_that("Parallel computation gives the same results as sequential", {
    skip_on_cran()

    if (getOption("dtwclust_skip_par_tests", FALSE))
        skip("Parallel tests disabled explicitly.")

    cat(" (with", num_workers, "workers)\n")

    require(doParallel)

    cl <- makeCluster(num_workers)
    clusterExport(cl, "num_workers", environment())
    invisible(clusterEvalQ(cl, {
        library(dtwclust)
        # environment variables get inherited by the workers when they are created, so set this
        if (!nzchar(Sys.getenv("R_COVR"))) RcppParallel::setThreadOptions(num_workers)
        # to test that other RNGkinds won't affect
        RNGkind("default")
    }))
    registerDoParallel(cl)

    # Filter excludes files that have "parallel" in them, otherwise it would be recursive
    res <- test_dir("./", filter = "parallel", invert = TRUE)

    sapply(clusterEvalQ(cl, RNGkind()[1L]), function(current_rngkind) {
        expect_identical(current_rngkind, default_rngkind)
    })

    stopCluster(cl)
    stopImplicitCluster()
    registerDoSEQ()
    rm(cl)
})

test_that("Parallel FORK computation gives the same results as sequential", {
    skip_on_cran()
    skip_on_travis()
    skip_on_os("windows")
    skip_if(nzchar(Sys.getenv("R_COVR")), "calculating coverage")

    if (getOption("dtwclust_skip_par_tests", FALSE))
        skip("Parallel tests disabled explicitly.")

    # Also test FORK in Linux
    cat(" - Test FORKs:\n")

    cl <- makeCluster(num_workers - 1L, "FORK")
    invisible(clusterEvalQ(cl, RcppParallel::setThreadOptions(2L)))
    registerDoParallel(cl)

    # Filter excludes files that have "parallel" in them, otherwise it would be recursive
    res <- test_dir("./", filter = "parallel", invert = TRUE)
    expect_s3_class(res, "testthat_results")

    stopCluster(cl)
    stopImplicitCluster()
    registerDoSEQ()
})
