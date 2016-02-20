context("Test parallel")

# =================================================================================================
# Compare to sequential
# - This breaks 'CMD check' within RStudio, it just crashes when it creates the parallel workers
# - It still works if I manually build and check the project in the command line
# =================================================================================================

test_that("Parallel computation gives the same results as sequential", {
     skip_on_cran()

     cat("\n")

     require(doParallel)

     if (.Platform$OS.type == "windows")
          cl <- makeCluster(detectCores())
     else
          cl <- makeCluster(detectCores(), "FORK")

     invisible(clusterEvalQ(cl, library(dtwclust)))

     registerDoParallel(cl)

     ## Filter excludes files that have "parallel" in them, otherwise it would be recursive
     test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

     stopCluster(cl)
     stopImplicitCluster()
     registerDoSEQ()
})
