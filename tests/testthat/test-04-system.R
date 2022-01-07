source("system/invalid-inputs.R", TRUE)
source("system/data-formats.R", TRUE)
source("system/preproc.R", TRUE)
source("system/fuzzy.R", TRUE)
source("system/hierarchical.R", TRUE)
source("system/partitional.R", TRUE)
source("system/comparisons.R", TRUE)

test_that("The RNGkind was not affected by dtwclust.", {
    expect_identical(RNGkind()[1L], default_rngkind)
})
