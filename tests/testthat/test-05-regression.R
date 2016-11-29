context("Regression tests")

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    source("regression/proxy.R")
    source("regression/dtwb.R")
    source("regression/family-distmat.R")
    source("regression/family-centroids.R")
    source("regression/custom-dist.R")
    source("regression/cvis.R")
    source("regression/methods.R")
    source("regression/clusterings.R")
}
