context("    Compare clusterings")

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

cfgs_mats <- compare_clusterings_configs(types = "h", k = 2L,
                                         preprocs = pdc_configs(
                                             "preproc",
                                             none = list()
                                         ),
                                         distances = pdc_configs(
                                             "distance",
                                             gak = list(window.size = 20L,
                                                        sigma = 100)
                                         ),
                                         centroids = pdc_configs(
                                             "centroid",
                                             DBA = list(window.size = 20L,
                                                        max.iter = 5L)
                                         )
)

# ==================================================================================================
# Compare clusterings
# ==================================================================================================

test_that("Compare clusterings works for the minimum set with all possibilities.", {
    expect_warning(errorpass_comp <- compare_clusterings(data_subset, c("p", "h", "f"),
                                                         configs = compare_clusterings_configs(k = 2L:3L),
                                                         seed = 932L, return.objects = TRUE,
                                                         .errorhandling = "pass"),
                   "names")

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
                                                  pick.clus = function(...) stop("NO!"),
                                                  lbls = labels_subset),
                   "pick.clus")
    expect_null(no_pick$pick)
    expect_true(!is.null(no_pick$scores))

    expect_warning(compare_clusterings(data_reinterpolated_subset, c("f"),
                                       configs = cfgs, seed = 392L,
                                       score.clus = bad_score_fun,
                                       lbls = labels_subset),
                   "scores.*not.*appended")

    type_score <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                      configs = cfgs, seed = 392L,
                                      score.clus = type_score_fun,
                                      lbls = labels_subset)

    type_score_objs <- compare_clusterings(data_reinterpolated_subset, c("f"),
                                           configs = cfgs, seed = 392L,
                                           return.objects = TRUE,
                                           score.clus = type_score_fun,
                                           lbls = labels_subset)

    expect_identical(no_pick$results, type_score$results)
    expect_identical(no_pick$results, type_score_objs$results)

    expect_output(
        all_comparisons <- compare_clusterings(data_reinterpolated_subset,
                                               c("p", "h", "f", "t"),
                                               configs = cfgs, seed = 392L,
                                               trace = TRUE,
                                               score.clus = score_fun,
                                               pick.clus = pick_fun,
                                               return.objects = TRUE,
                                               shuffle.configs = TRUE,
                                               lbls = labels_subset)
    )

    gak_comparison <- compare_clusterings(data_subset, "p",
                                          configs = cfgs_gak, seed = 190L,
                                          score.clus = score_fun,
                                          lbls = labels_subset)

    do_this <- if (foreach::getDoParWorkers() > 1L) base::eval else testthat::expect_warning
    do_this({
        dba_comparison <- compare_clusterings(data_multivariate, "h",
                                              configs = cfgs_dba, seed = 294L,
                                              score.clus = score_fun,
                                              lbls = labels_subset)
    })

    sdtwc_comparison <- compare_clusterings(data_subset, "h",
                                            configs = cfgs_sdtwc, seed = 3290L,
                                            score.clus = score_fun,
                                            lbls = labels_subset)

    ## rds
    all_comparisons$pick <- reset_nondeterministic(all_comparisons$pick)
    all_comparisons$pick@call <- call("zas", foo = "bar")
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
})

# ==================================================================================================
# clean
# ==================================================================================================

rm(list = setdiff(ls(), ols))
