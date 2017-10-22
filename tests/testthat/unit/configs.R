context("\tConfigs")

# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

# ==================================================================================================
# pdc_configs
# ==================================================================================================

test_that("Using pdc_configs() works as expected.", {
    expect_warning(pdc_configs("d", tadpole = list(a = "a")))

    share_base <- c("partitional", "hierarchical", "fuzzy", "tadpole")

    pdc_p <- list(pp = list(p1 = 1L))
    pdc_h <- list(ph = list(h1 = 1L, h2 = 1L:2L))
    pdc_f <- list(pf = list(f1 = "a", f2 = c(FALSE, TRUE), f3 = 1:3))
    pdc_t <- list(pt = list(t2 = c("x", "y")))

    pdc_nrows <- c(
        partitional = 2L,
        hierarchical = 3L,
        fuzzy = 7L,
        tadpole = 3L
    )

    pdc_args <- list(
        # shared
        none = list(),
        # specific
        partitional = pdc_p,
        hierarchical = pdc_h,
        fuzzy = pdc_f,
        tadpole = pdc_t
    )

    nrows_is_correct <- times(100L) %do% {
        this_args <- sample(pdc_args)[1L:sample(length(pdc_args), 1L)]
        this_share <- sample(share_base)[1L:sample(length(share_base), 1L)]
        type <- sample(c("preproc", "distance", "centroid"), 1L)

        if (type == "distance") {
            this_args$tadpole <- NULL
            this_share <- this_share[setdiff(names(this_share), "tadpole")]
        }

        if (length(this_args) == 0L || length(this_share) == 0L)
            return(TRUE)

        this_cfg <- do.call(pdc_configs, c(
            this_args,
            list(
                type = type,
                share.config = this_share
            )))

        decrease <- c(
            partitional = 0L,
            hierarchical = 0L,
            fuzzy = 0L,
            tadpole = 0L
        )

        if (is.null(this_args$none)) {
            decrease <- decrease + 1L

        } else {
            must_decrease <- setdiff(share_base, this_share)
            decrease[must_decrease] <- decrease[must_decrease] + 1L
        }

        must_decrease <- sapply(this_args[names(decrease)], is.null)
        decrease[must_decrease] <- decrease[must_decrease] + pdc_nrows[must_decrease] - 1L
        expected_nrows <- pdc_nrows - decrease
        all(expected_nrows[names(this_cfg)] == sapply(this_cfg, nrow))
    }

    expect_true(all(nrows_is_correct))
})

# ==================================================================================================
# compare_clusterings_configs
# ==================================================================================================

test_that("Using compare_clusterings_configs() works as expected.", {
    ## ---------------------------------------------------------- controls errors
    expect_error(compare_clusterings_configs(controls = "a"),
                 "controls.*list")
    expect_error(compare_clusterings_configs(controls = list(partitional_control())),
                 "controls.*name")
    expect_error(compare_clusterings_configs(controls = list(partitional = partitional_control())),
                 "controls.*types")

    ## ---------------------------------------------------------- preprocs errors
    expect_error(compare_clusterings_configs(preprocs = "a"),
                 "preprocs.*list")
    expect_error(compare_clusterings_configs(preprocs = list(pdc_configs())),
                 "preprocs.*name")
    expect_error(compare_clusterings_configs(preprocs = list(partitional = pdc_configs())),
                 "preprocs.*types")

    ## ---------------------------------------------------------- distances errors
    expect_error(compare_clusterings_configs(distances = "a"),
                 "distances.*list")
    expect_error(compare_clusterings_configs(distances = list(pdc_configs())),
                 "distances.*name")
    expect_error(compare_clusterings_configs(distances = list(partitional = pdc_configs())),
                 "distances.*types")

    ## ---------------------------------------------------------- centroids errors
    expect_error(compare_clusterings_configs(centroids = "a"),
                 "centroids.*list")
    expect_error(compare_clusterings_configs(centroids = list(pdc_configs())),
                 "centroids.*name")
    expect_error(compare_clusterings_configs(centroids = list(partitional = pdc_configs())),
                 "centroids.*types")

    ## ---------------------------------------------------------- working configs
    nrows_each_config <- sapply(compare_clusterings_configs(), nrow)
    expect_true(all(nrows_each_config == 1L))

    nrows_each_custom_config <- sapply(FUN = nrow, X = compare_clusterings_configs(
        types = c("p", "h", "f", "t"),
        controls = list(
            partitional = partitional_control(),
            hierarchical = hierarchical_control(),
            fuzzy = fuzzy_control(),
            tadpole = tadpole_control(dc = 1.5, window.size = 1L)
        ),
        preprocs = pdc_configs("p", none = list(foo = "bar")),
        distances = pdc_configs("d", sbd = list(foo = "bar")),
        centroids = pdc_configs("c", default = list(foo = "bar"))
    ))
    expect_true(all(nrows_each_custom_config == 1L))
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
