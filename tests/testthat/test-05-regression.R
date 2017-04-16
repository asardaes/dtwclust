context("Regression tests")

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    source("regression/proxy.R", TRUE)
    source("regression/dtwb.R", TRUE)
    source("regression/family-distmat.R", TRUE)
    source("regression/family-centroids.R", TRUE)
    source("regression/custom-dist.R", TRUE)
    source("regression/cvis.R", TRUE)
    source("regression/methods.R", TRUE)
    source("regression/clusterings.R", TRUE)
    source("regression/comparisons.R", TRUE)
}
