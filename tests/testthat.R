library(dtwclust)
library(testthat)

#' To test in a local machine:
#' Sys.setenv(NOT_CRAN = "true"); test_dir("tests/testthat/")
#' OR
#' devtools::test()
#'
#' To disable parallel tests, before calling test() run:
#'
#' options(skip_par_tests = TRUE)
#'
testthat::test_check("dtwclust")
