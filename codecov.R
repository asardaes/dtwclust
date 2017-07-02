library(covr)
codecov(type = "all",
        line_exclusions = list(
            "R/pkg.R",
            "R/create-dtwclust.R",
            "R/dtwclust-classes.R",
            "R/dtwclust-methods.R",
            "R/compute-envelope.R" = 55L:58L,
            "R/partitional.R" = 130L:253L
        ))
