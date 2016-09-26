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

     cat("\n")

     require(doParallel)

     cl <- makeCluster(detectCores())
     invisible(clusterEvalQ(cl, library(dtwclust)))
     registerDoParallel(cl)

     ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
     test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

     stopCluster(cl)
     stopImplicitCluster()
     registerDoSEQ()

     ## Also test FORK in Unix
     if (.Platform$OS.type != "windows") {
          cat("Test FORKs:\n")

          rm(cl)
          cl <- makeCluster(detectCores(), "FORK")
          registerDoParallel(cl)

          ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
          test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

          stopCluster(cl)
          stopImplicitCluster()
          registerDoSEQ()
     }
})
