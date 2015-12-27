context("Test parallel")

# =================================================================================================
# Compare to sequential
# =================================================================================================

test_that("Parallel computation gives the same results as sequential", {
     skip_on_cran()

     require(doParallel)

     if (.Platform$OS.type == "windows") {
          cl <- makeCluster(detectCores())
          invisible(clusterEvalQ(cl, library(dtwclust)))

     } else {
          cl <- makeCluster(detectCores(), "FORK")
     }

     registerDoParallel(cl)

     ## The 'distmat' slot in a couple of experiments changes because the parallel version forces
     ## symmetry by only calculating the upper triangular...
     test_dir("./", filter = "^(?!.*parallel).*$", perl = TRUE)

     stopCluster(cl)
     stopImplicitCluster()
     registerDoSEQ()
})
