## Fuzzy preprocessing: calculate autocorrelation up to 50th lag
acf_fun <- function(dat, ...) {
    lapply(dat, function(x) {
        as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)
    })
}

## Define overall configuration
cfgs <- compare_clusterings_configs(
    types = c("p", "h", "f", "t"),
    k = 19L:20L,
    controls = list(
        partitional = partitional_control(
            iter.max = 30L,
            nrep = 1L
        ),
        hierarchical = hierarchical_control(
            method = "all"
        ),
        fuzzy = fuzzy_control(
            ## notice the vector
            fuzziness = c(2, 2.5),
            iter.max = 30L
        ),
        tadpole = tadpole_control(
            ## notice the vectors
            dc = c(1.5, 2),
            window.size = 19L:20L
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
        sbd = list(),
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
            default = list(),
            shape_extraction = list(znorm = TRUE)
        )
    )
)

## Number of configurations is returned as attribute
num_configs <- sapply(cfgs, attr, which = "num.configs")
cat("\nTotal number of configurations without considering optimizations:",
    sum(num_configs),
    "\n\n")

## Define evaluation function based on CVI: Variation of Information
score_fun <- function(obj_list, ...) {
    sapply(obj_list, function(obj) {
        cvi(obj@cluster, CharTrajLabels, type = "VI")
    })
}

## Function that chooses best result
pick_fun <- function(scores, obj_lists, ...) {
    best_considering_type <- sapply(scores, which.min)
    best_overall <- which.min(mapply(scores, best_considering_type,
                                     FUN = function(score, id) { score[id] }))

    best_obj <- obj_lists[[best_overall]][[best_considering_type[best_overall]]]

    ## return
    best_obj
}

# ====================================================================================
# Short run with only fuzzy clustering
# ====================================================================================

comparison_short <- compare_clusterings(CharTraj, types = c("f"), configs = cfgs,
                                        seed = 293L, trace = TRUE,
                                        score.clus = score_fun, pick.clus = pick_fun,
                                        return.objects = TRUE)

\dontrun{
# ====================================================================================
# Parallel run with all comparisons
# ====================================================================================

require(doParallel)
registerDoParallel(cl <- makeCluster(detectCores()))

comparison_long <- compare_clusterings(CharTraj, types = c("p", "h", "f", "t"),
                                       configs = cfgs,
                                       seed = 293L, trace = TRUE,
                                       score.clus = score_fun,
                                       pick.clus = pick_fun,
                                       return.objects = TRUE)

# ------------------------------------------------------------------------------------
# Using all external CVIs and majority vote
# ------------------------------------------------------------------------------------

score_external <- function(obj_list, ...) {
    scores <- lapply(obj_list, function(obj) {
        indices <- cvi(obj@cluster, CharTrajLabels, type = "external")

        ## invert VI to consider maximization
        indices["VI"] <- 1 / indices["VI"]

        ## return
        indices
    })

    ## return
    do.call(rbind, scores)
}

pick_majority <- function(scores, obj_lists, ...) {
    majority <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }

    best_considering_type <- sapply(scores, function(score) {
        best_by_cvi <- apply(score, 2L, which.max)
        majority(best_by_cvi)
    })

    best_cvis_overall <- do.call(rbind,
                                 mapply(scores, best_considering_type,
                                        SIMPLIFY = FALSE,
                                        FUN = function(score, row_id) {
                                            score[row_id, , drop = FALSE]
                                        }))

    best_overall <- majority(apply(best_cvis_overall, 2L, which.max))

    best_obj <- obj_lists[[best_overall]][[best_considering_type[best_overall]]]

    ## to find config later, see 'best config' below
    attr(best_obj, "config_id") <- c(best_overall,
                                     best_considering_type[best_overall])

    ## return
    best_obj
}

comparison_majority <- compare_clusterings(CharTraj, types = c("p", "h", "f", "t"),
                                           configs = cfgs,
                                           seed = 84L, trace = TRUE,
                                           score.clus = score_external,
                                           pick.clus = pick_majority,
                                           return.objects = TRUE)

plot(comparison_majority$pick)

## best config
config_id <- attr(comparison_majority$pick, "config_id")
print(comparison_majority$results[[config_id[1L]]][config_id[2L], , drop = FALSE])

stopCluster(cl); registerDoSEQ()

# ====================================================================================
# A run with only partitional clusterings
# ====================================================================================

p_cfgs <- compare_clusterings_configs(types = "p", k = 19L:21L,
                                      controls = list(
                                          partitional = partitional_control(
                                              iter.max = 20L,
                                              nrep = 8L
                                          )
                                      ),
                                      preprocs = pdc_configs(
                                          "preproc",
                                          none = list(),
                                          zscore = list(center = c(FALSE, TRUE))
                                      ),
                                      distances = pdc_configs(
                                          "distance",
                                          sbd = list(),
                                          dtw_basic = list(window.size = 19L:20L,
                                                           norm = c("L1", "L2")),
                                          gak = list(window.size = 19L:20L,
                                                     sigma = 100)
                                      ),
                                      centroids = pdc_configs(
                                          "centroid",
                                          partitional = list(
                                              pam = list(),
                                              shape = list()
                                          )
                                      )
)

# Remove redundant (shape centroid always uses zscore preprocessing)
id_redundant <- p_cfgs$partitional$preproc == "none" &
    p_cfgs$partitional$centroid == "shape"
p_cfgs$partitional <- p_cfgs$partitional[!id_redundant, ]

# LONG! 20 minutes or so, sequentially
comparison_partitional <- compare_clusterings(CharTraj, types = "p",
                                              configs = p_cfgs,
                                              seed = 32903L, trace = TRUE,
                                              score.clus = score_fun,
                                              pick.clus = pick_fun,
                                              shuffle.configs = TRUE,
                                              return.objects = TRUE)
}
