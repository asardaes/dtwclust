library(dtwclust)
library(foreach)
library(testthat)

# coverage for multi-threading might not be possible (?)
if (nzchar(Sys.getenv("R_COVR"))) RcppParallel::setThreadOptions(1L)
# old reporter for CMD checks
options(testthat.default_reporter = "summary")

#' To test in a local machine:
#' Sys.setenv(NOT_CRAN = "true"); test_dir("tests/testthat/")
#' OR
#' devtools::test() # broken, can't figure out why
#'
#' To disable parallel tests, before calling test() run:
#'
#'   options(dtwclust_skip_par_tests = TRUE)
#'
testthat::test_check("dtwclust")
