context("\tCompare clusterings")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

acf_fun <- function(dat, ...) {
    lapply(dat, function(x) {
        as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)
    })
}

score_fun <- function(obj_list, lbls, ...) {
    sapply(obj_list, function(obj) {
        cvi(obj@cluster, lbls, type = "VI")
    })
}

pick_fun <- function(scores, obj_lists, ...) {
    best_considering_type <- sapply(scores, which.min)
    best_overall <- which.min(mapply(scores, best_considering_type,
                                     FUN = function(score, id) { score[id] }))

    best_obj <- obj_lists[[best_overall]][[best_considering_type[best_overall]]]

    ## return
    best_obj
}

type_score_fun <- list(fuzzy = function(obj_list, lbls, ...) {
    sapply(obj_list, function(obj) {
        cvi(obj@cluster, lbls, type = "VI")
    })
})

cfgs <- compare_clusterings_configs(c("p", "h", "f", "t"), k = 2L:3L,
                                    controls = list(
                                        partitional = partitional_control(
                                            pam.precompute = c(FALSE, TRUE),
                                            iter.max = 10L,
                                            nrep = 2L
                                        ),
                                        hierarchical = hierarchical_control(
                                            method = "all"
                                        ),
                                        fuzzy = fuzzy_control(
                                            fuzziness = c(2, 2.5),
                                            iter.max = 10L,
                                            delta = c(0.1, 0.01)
                                        ),
                                        tadpole = tadpole_control(
                                            dc = c(1.5, 2),
                                            window.size = 19L:20L,
                                            lb = c("lbk", "lbi")
                                        )
                                    ),
                                    preprocs = pdc_configs(
                                        type = "preproc",
                                        ## shared
                                        none = list(),
                                        zscore = list(center = c(FALSE)),
                                        ## only for fuzzy
                                        fuzzy = list(
                                            acf_fun = list()
                                        ),
                                        ## only for tadpole
                                        tadpole = list(
                                            reinterpolate = list(new.length = 205L)
                                        ),
                                        ## specify which should consider the shared ones
                                        share.config = c("p", "h")
                                    ),
                                    distances = pdc_configs(
                                        type = "distance",
                                        dtw_basic = list(
                                            norm = c("L1", "L2"),
                                            window.size = 18L
                                        ),
                                        fuzzy = list(
                                            L2 = list()
                                        ),
                                        share.config = c("p", "h")
                                    ),
                                    centroids = pdc_configs(
                                        type = "centroid",
                                        partitional = list(
                                            pam = list()
                                        ),
                                        ## special name 'default'
                                        hierarchical = list(
                                            default = list()
                                        ),
                                        fuzzy = list(
                                            fcmdd = list()
                                        ),
                                        tadpole = list(
                                            shape_extraction = list(znorm = TRUE)
                                        )
                                    )
)

# =================================================================================================
# Compare clusterings
# =================================================================================================

test_that("Compare clusterings works for the minimum set with all possibilities.", {
    expect_warning(no_score <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                                   configs = cfgs, seed = 392L))
    expect_null(no_score$scores)

    expect_warning(no_pick <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                                  configs = cfgs, seed = 392L,
                                                  score.clus = score_fun,
                                                  lbls = labels_subset))
    expect_null(no_pick$pick)
    expect_true(!is.null(no_pick$scores))

    expect_warning(type_score <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                                     configs = cfgs, seed = 392L,
                                                     score.clus = type_score_fun,
                                                     lbls = labels_subset))

    expect_identical(no_pick$results, type_score$results)

    mute <- capture.output(all_comparisons <- compare_clusterings(data_reinterpolated_subset,
                                                                  c("p", "h", "f", "t"),
                                                                  configs = cfgs, seed = 392L,
                                                                  trace = TRUE,
                                                                  score.clus = score_fun,
                                                                  pick.clus = pick_fun,
                                                                  lbls = labels_subset))

    ## rds
    all_comparisons$pick <- reset_nondeterministic(all_comparisons$pick)
    all_comparisons$pick@call <- call("zas", foo = "bar")
    all_comparisons$proc_time <- NULL

    assign("all_comp", all_comparisons, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
