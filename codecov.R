library(covr)
options(covr.exclude_pattern = TRUE, covr.exclude_start = TRUE, covr.exclude_end = TRUE)
codecov(type = "all",
        line_exclusions = list(
            "R/pkg.R"
        ))
