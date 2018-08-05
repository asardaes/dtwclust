library("covr")
covr_flags <- getOption("covr.flags")
options(covr.exclude_pattern = rex::rex("#" %or% "//", any_spaces, "nocov"),
        covr.exclude_start = rex::rex("#" %or% "//", any_spaces, "nocov", any_spaces, "start"),
        covr.exclude_end = rex::rex("#" %or% "//", any_spaces, "nocov", any_spaces, "end"),
        covr.gcov = "gcov-7", # see https://github.com/r-lib/covr/issues/176
        covr.flags = sapply(covr_flags, function(dummy) { "--coverage" }))
covr::codecov(type = "tests", quiet = FALSE)

#' to run locally, comment the line with covr.gcov above, run all code above until codecov(), and:
#'
#' Sys.setenv(NOT_CRAN = "true"); covr::report(package_coverage(type = "tests", quiet = FALSE))
