library(covr)
options(covr.exclude_pattern = TRUE, covr.exclude_start = TRUE, covr.exclude_end = TRUE)
codecov(type = "all",
        line_exclusions = list(
            "R/pkg.R",
            "R/tslist.R",
            "R/dtwclust.R",
            "R/all-cent.R",
            "R/ddist.R",
            "R/create-dtwclust.R",
            "R/dtwclust-classes.R",
            "R/dtwclust-methods.R",
            "R/compute-envelope.R" = 55L:58L,
            "R/partitional.R" = 130L:253L
        ))
