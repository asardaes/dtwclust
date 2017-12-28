library("covr")
options(covr.exclude_pattern = rex::rex("#" %or% "//", any_spaces, "nocov"),
        covr.exclude_start = rex::rex("#" %or% "//", any_spaces, "nocov", any_spaces, "start"),
        covr.exclude_end = rex::rex("#" %or% "//", any_spaces, "nocov", any_spaces, "end"))
covr::codecov(type = "tests", line_exclusions = list(
    "R/pkg.R",
    "src/dtwclust.h",
    "src/dtwclust++.h"
))
