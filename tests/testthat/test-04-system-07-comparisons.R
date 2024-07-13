# ==================================================================================================
# setup
# ==================================================================================================

## Original objects in env
ols <- ls()

acf_fun <- function(dat, ...) {
    lapply(dat, function(x) {
        as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)
    })
}

evaluators <- suppressMessages(cvi_evaluators("VI", ground.truth = labels_subset))
score_fun <- evaluators$score
pick_fun <- evaluators$pick
type_score_fun <- list(fuzzy = score_fun)
bad_score_fun <- list(fuzzy = function(...) { list(1L:2L, 3L:5L) })

cfgs <- compare_clusterings_configs(c("p", "h", "f", "t"), k = 2L:3L,
                                    controls = list(
                                        partitional = partitional_control(
                                            pam.precompute = c(FALSE, TRUE),
                                            iter.max = 10L,
                                            nrep = 2L,
                                            version = 2L
                                        ),
                                        hierarchical = hierarchical_control(
                                            method = "all"
                                        ),
                                        fuzzy = fuzzy_control(
                                            fuzziness = c(2, 2.5),
                                            iter.max = 10L,
                                            delta = c(0.1, 0.01),
                                            version = 2L
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

cfgs_gak <- compare_clusterings_configs(types = "p", k = 2L:3L,
                                        controls = list(
                                            partitional = partitional_control(
                                                iter.max = 5L,
                                                nrep = 1L,
                                                version = 2L
                                            )
                                        ),
                                        preprocs = pdc_configs(
                                            "preproc",
                                            none = list()
                                        ),
                                        distances = pdc_configs(
                                            "distance",
                                            gak = list(window.size = 20L, sigma = c(100, 120))
                                        ),
                                        centroids = pdc_configs(
                                            "centroid",
                                            pam = list()
                                        )
)

cfgs_dba <- compare_clusterings_configs(types = "h", k = 2L:3L,
                                        preprocs = pdc_configs(
                                            "preproc",
                                            none = list()
                                        ),
                                        distances = pdc_configs(
                                            "distance",
                                            dtw_basic = list(window.size = 20L)
                                        ),
                                        centroids = pdc_configs(
                                            "centroid",
                                            DBA = list(window.size = 20L,
                                                       max.iter = 5L)
                                        )
)

cfgs_sdtwc <- compare_clusterings_configs(types = "h", k = 2L,
                                          preprocs = pdc_configs(
                                              "preproc",
                                              none = list()
                                          ),
                                          distances = pdc_configs(
                                              "distance",
                                              sdtw = list()
                                          ),
                                          centroids = pdc_configs(
                                              "centroid",
                                              sdtw_cent = list()
                                          )
)

# ==================================================================================================
# CVI evaluators
# ==================================================================================================

test_that("cvi_evaluators work nicely with compare_clusterings.", {
    # these should be unit tests, but oh well
    expect_error(cvi_evaluators("external"), "ground.truth")
    expect_error(cvi_evaluators("VI"), "ground.truth")

    pick_inconclusive <- cvi_evaluators(c("RI", "ARI"), ground.truth = labels_subset)$pick

    res <- list(data.frame(RI = 0:1, ARI = 1:0))
    expect_error(pick_inconclusive(res), "inconclusive")

    res <- list(data.frame(RI = c(0,1), ARI = c(0,2)),
                data.frame(RI = c(2,0), ARI = c(1,0)))
    expect_error(pick_inconclusive(res), "inconclusive")

    # smaller cfgs
    cfgs <- compare_clusterings_configs(c("h", "f"), k = 2L:3L,
                                        controls = list(
                                            hierarchical = hierarchical_control(
                                                method = "average"
                                            ),
                                            fuzzy = fuzzy_control(
                                                fuzziness = c(2, 2.5),
                                                iter.max = 10L
                                            )
                                        ),
                                        preprocs = pdc_configs(
                                            type = "preproc",
                                            ## shared
                                            none = list(),
                                            ## only for fuzzy
                                            fuzzy = list(
                                                acf_fun = list()
                                            ),
                                            ## specify which should consider the shared ones
                                            share.config = c("h")
                                        ),
                                        distances = pdc_configs(
                                            type = "distance",
                                            sbd = list(),
                                            fuzzy = list(
                                                L2 = list()
                                            ),
                                            share.config = c("h")
                                        ),
                                        centroids = pdc_configs(
                                            type = "centroid",
                                            ## special name 'default'
                                            hierarchical = list(
                                                default = list()
                                            ),
                                            fuzzy = list(
                                                fcmdd = list()
                                            )
                                        )
    )

    # non-fuzzy
    suppressMessages({
        ev_valid <- cvi_evaluators("valid", ground.truth = labels_subset)
        ev_internal <- cvi_evaluators("internal", ground.truth = labels_subset)
        ev_external <- cvi_evaluators("external", ground.truth = labels_subset)
        ev_vi <- cvi_evaluators("VI", ground.truth = labels_subset)
    })

    # valid
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("h"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_valid$score,
                                     pick.clus = ev_valid$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("h"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_valid$score,
                                        pick.clus = ev_valid$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # internal
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("h"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_internal$score,
                                     pick.clus = ev_internal$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("h"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_internal$score,
                                        pick.clus = ev_internal$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # external
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("h"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_external$score,
                                     pick.clus = ev_external$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("h"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_external$score,
                                        pick.clus = ev_external$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # vi
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("h"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_vi$score,
                                     pick.clus = ev_vi$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("h"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_vi$score,
                                        pick.clus = ev_vi$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # fuzzy
    suppressMessages({
        ev_valid <- cvi_evaluators("valid", TRUE, ground.truth = labels_subset)
        ev_internal <- cvi_evaluators("internal", TRUE, ground.truth = labels_subset)
        ev_external <- cvi_evaluators("external", TRUE, ground.truth = labels_subset)
        ev_vi <- cvi_evaluators("VI", TRUE, ground.truth = labels_subset)
    })

    # valid
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("f"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_valid$score,
                                     pick.clus = ev_valid$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("f"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_valid$score,
                                        pick.clus = ev_valid$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # internal
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("f"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_internal$score,
                                     pick.clus = ev_internal$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("f"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_internal$score,
                                        pick.clus = ev_internal$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # external
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("f"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_external$score,
                                     pick.clus = ev_external$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("f"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_external$score,
                                        pick.clus = ev_external$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)

    # vi
    with_objs <- compare_clusterings(data_reinterpolated_subset,
                                     c("f"),
                                     configs = cfgs, seed = 392L,
                                     score.clus = ev_vi$score,
                                     pick.clus = ev_vi$pick,
                                     return.objects = TRUE)
    without_objs <- compare_clusterings(data_reinterpolated_subset,
                                        c("f"),
                                        configs = cfgs, seed = 392L,
                                        score.clus = ev_vi$score,
                                        pick.clus = ev_vi$pick,
                                        return.objects = FALSE)
    expect_identical(with_objs$results, without_objs$results)
    expect_identical(with_objs$pick$config, without_objs$pick)
})

# ==================================================================================================
# Compare clusterings
# ==================================================================================================

test_that("Compare clusterings works for the minimum set with all possibilities.", {
    errored_cfg <- compare_clusterings_configs("t",controls = list(
        tadpole = tadpole_control(1.5, 10L)
    ))
    errored_cfg$tadpole$window.size <- NULL
    expect_warning(
        expect_warning(
            errored <- compare_clusterings(data_reinterpolated_subset, "t", errored_cfg,
                                           .errorhandling = "pass",
                                           score.clus = function(...) {})
        )
    )
    expect_true(inherits(errored$scores$tadpole[[1L]], "error"))

    expect_warning(
        expect_warning(
            errorpass_comp <- compare_clusterings(data_subset, c("p", "h", "f"),
                                                  configs = compare_clusterings_configs(k = 2L:3L),
                                                  seed = 932L, return.objects = TRUE,
                                                  .errorhandling = "pass"),
            "names"
        )
    )

    expect_true(inherits(errorpass_comp$objects.fuzzy[[1L]], "error"))

    expect_warning(errorrm_comp <- compare_clusterings(data_subset, c("p", "h", "f"),
                                                       configs = compare_clusterings_configs(),
                                                       seed = 932L, return.objects = TRUE,
                                                       .errorhandling = "remove"))

    expect_null(errorrm_comp$objects.fuzzy)

    expect_warning(no_score <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                                   configs = cfgs, seed = 392L,
                                                   return.objects = TRUE,
                                                   score.clus = function(...) stop("NO!")),
                   "score.clus")
    expect_null(no_score$scores)

    expect_warning(no_pick <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                                  configs = cfgs, seed = 392L,
                                                  score.clus = score_fun,
                                                  pick.clus = function(...) stop("NO!")),
                   "pick.clus")
    expect_null(no_pick$pick)
    expect_true(!is.null(no_pick$scores$fuzzy))

    expect_warning(no_pick_with_objects <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                                               configs = cfgs, seed = 392L,
                                                               return.objects = TRUE,
                                                               score.clus = score_fun,
                                                               pick.clus = function(...) stop("NO!")),
                   "pick.clus")
    expect_null(no_pick_with_objects$pick)
    expect_true(!is.null(no_pick_with_objects$scores$fuzzy))

    expect_warning(compare_clusterings(data_reinterpolated_subset, c("f"),
                                       configs = cfgs, seed = 392L,
                                       score.clus = bad_score_fun),
                   "scores.*not.*appended")

    type_score <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                      configs = cfgs, seed = 392L,
                                      score.clus = type_score_fun)

    type_score_objs <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                           configs = cfgs, seed = 392L,
                                           return.objects = TRUE,
                                           score.clus = type_score_fun)

    expect_identical(no_pick$results, type_score$results)
    expect_identical(no_pick$results, type_score_objs$results)

    expect_output(
        suppressMessages(
            all_comparisons <- compare_clusterings(data_reinterpolated_subset,
                                                   c("p", "h", "f", "t"),
                                                   configs = cfgs, seed = 392L,
                                                   trace = TRUE,
                                                   score.clus = score_fun,
                                                   pick.clus = pick_fun,
                                                   return.objects = TRUE,
                                                   shuffle.configs = TRUE)
        )
    )

    expect_equal_slots(
        repeat_clustering(data_reinterpolated_subset, all_comparisons, "config3_1"),
        all_comparisons$objects.partitional$config3_1
    )

    expect_equal_slots(
        repeat_clustering(data_reinterpolated_subset, all_comparisons, "config10_1"),
        all_comparisons$objects.hierarchical$config10_1
    )

    expect_equal_slots(
        repeat_clustering(data_reinterpolated_subset, all_comparisons, "config18_2"),
        all_comparisons$objects.tadpole$config18_2
    )

    gak_comparison <- compare_clusterings(data_subset, "p",
                                          configs = cfgs_gak, seed = 190L,
                                          score.clus = score_fun)

    do_this <- if (foreach::getDoParWorkers() > 1L) base::eval else testthat::expect_warning
    do_this({
        dba_comparison <- compare_clusterings(data_multivariate, "h",
                                              configs = cfgs_dba, seed = 294L,
                                              score.clus = score_fun)
    })

    suppressWarnings(
        sdtwc_comparison <- compare_clusterings(data_subset, "h",
                                                configs = cfgs_sdtwc, seed = 3290L,
                                                score.clus = score_fun)
    )

    ## rds
    all_comparisons$pick$object <- reset_nondeterministic(all_comparisons$pick$object)
    all_comparisons$pick$object@call <- call("zas", foo = "bar")
    all_comparisons$proc_time <- NULL
    all_comparisons$objects.partitional <- NULL
    all_comparisons$objects.hierarchical <- NULL
    all_comparisons$objects.fuzzy <- NULL
    all_comparisons$objects.tadpole <- NULL
    gak_comparison$proc_time <- NULL
    dba_comparison$proc_time <- NULL
    sdtwc_comparison$proc_time <- NULL

    assign("comp_all", all_comparisons, persistent)
    assign("comp_gak", gak_comparison, persistent)
    assign("comp_dba", dba_comparison, persistent)
    assign("comp_sdtwc", sdtwc_comparison, persistent)

    simple_gak <- cfgs_gak
    simple_gak$partitional$k <- 2L
    simple_gak <- compare_clusterings(data_subset, "p", simple_gak,
                                      seed = 3289L, score.clus = score_fun)
    expect_error(repeat_clustering(data_subset, simple_gak, "config1_1"))
    expect_s4_class(repeat_clustering(data_subset, simple_gak, "config1"), "TSClusters")
})

test_that("Results data frame for hierarchical clustering is correct (GH issue #57).", {
    cfgs <- compare_clusterings_configs(types = "h", k = 2L:3L,
                                        controls = list(
                                            hierarchical = hierarchical_control(method = "all")
                                        ),
                                        preprocs = pdc_configs(
                                            "preproc",
                                            none = list()
                                        ),
                                        distances = pdc_configs(
                                            "distance",
                                            sbd = list()
                                        ),
                                        centroids = pdc_configs(
                                            "centroid",
                                            shape_extraction = list()
                                        ))

    cmp <- compare_clusterings(data_subset, "h", configs = cfgs, seed = 329L, return.objects = TRUE)

    ignored <- Map(seq_len(nrow(cmp$results$hierarchical)),
                   split.data.frame(cmp$results$hierarchical, factor(cmp$results$hierarchical$config_id,
                                                                     cmp$results$hierarchical$config_id)), # prevent factor re-ordering
                   cmp$objects.hierarchical,
                   f = function(i, res, obj) {
                       expect_identical(obj@k, res$k, info = paste("Row", i))
                       expect_identical(obj@method, res$method, info = paste("Row", i))
                   })
})

test_that("step.pattern can be provided in a nested list (GH issue #59)", {
    cfgs <- compare_clusterings_configs(
        types = "p",
        k = 4L,
        controls = list(
            partitional = partitional_control(
                iter.max = 30L,
                nrep = 1L
            )
        ),
        preprocs = pdc_configs(
            type = "preproc",
            none = list(),
            share.config = c("p")
        ),
        distances = pdc_configs(
            type = "distance",
            dtw_basic = list(
                norm = c("L1", "L2"),
                step.pattern = list(dtw::symmetric1, dtw::symmetric2)
            ),
            share.config = c("p", "h")
        ),
        centroids = pdc_configs(
            type = "centroid",
            partitional = list(
                pam = list()
            )
        )
    )

    cmp <- compare_clusterings(data_subset, "p", configs = cfgs, seed = 329L, return.objects = TRUE)
    expect_identical(length(cmp$objects.partitional), 4L)
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
