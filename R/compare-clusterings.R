#' Helper function for [compare_clusterings_configs()]
#'
#' Create preprocessing, distance and centroid configurations for [compare_clusterings_configs()].
#' All functions use [base::expand.grid()].
#'
#' @export
#'
#' @param type Which type of function is being targeted by this configuration.
#' @param ... Any number of named nested lists with functions and arguments that will be shared by
#'   all clusterings. See details.
#' @param partitional A named list of lists with functions and arguments for partitional
#'   clusterings.
#' @param hierarchical A named list of lists with functions and arguments for hierarchical
#'   clusterings.
#' @param fuzzy A named list of lists with functions and arguments for fuzzy clusterings.
#' @param tadpole A named list of lists with functions and arguments for TADPole clusterings.
#' @param share.config A character vector specifying which clusterings should include the shared
#'   lists (the ones specified in `...`). It must be any combination of (possibly abbreviated):
#'   partitional, hierarchical, fuzzy, tadpole.
#'
#' @details
#'
#' The named lists are interpreted in the following way: the name of each element of the list will
#' be considered to be a function name, and the elements of the nested list will be the possible
#' parameters for the function. Each function must have at least an empty list. The parameters may
#' be vectors that specify different values to be tested.
#'
#' For preprocessing, the special name `none` signifies no preprocessing.
#'
#' For centroids, the special name `default` leaves the centroid unspecified.
#'
#' Please see the examples in [compare_clusterings()] to see how this is used.
#'
#' @return
#'
#' A list for each clustering, each of which includes a data frame with the computed configurations.
#'
pdc_configs <- function(type = c("preproc", "distance", "centroid"), ...,
                        partitional = NULL, hierarchical = NULL, fuzzy = NULL, tadpole = NULL,
                        share.config = c("p", "h", "f", "t"))
{
    type <- match.arg(type)

    shared <- list(...)

    specific <- list(partitional = partitional,
                     hierarchical = hierarchical,
                     fuzzy = fuzzy,
                     tadpole = tadpole)

    specific <- specific[!sapply(specific, is.null)]

    if (type == "distance") {
        if (!is.null(specific$tadpole)) warning("TADPole ignores distance configurations.")
        specific$tadpole <- NULL
    }

    share_missing <- missing(share.config)
    share.config <- match.arg(share.config, supported_clusterings, TRUE)

    if (type == "distance") {
        if (!share_missing && "tadpole" %in% share.config)
            warning("TADPole ignores distance configurations.")

        share.config <- setdiff(share.config, "tadpole")
    }

    ## =============================================================================================
    ## Shared configs
    ## =============================================================================================

    if (length(shared) > 0L && length(share.config) > 0L) {
        shared_names <- names(shared)

        shared_cfg <- mapply(shared, shared_names, SIMPLIFY = FALSE,
                             FUN = function(shared_args, fun) {
                                 cfg <- do.call(expand.grid,
                                                enlist(foo = fun,
                                                       dots = shared_args,
                                                       stringsAsFactors = FALSE))

                                 names(cfg)[1L] <- type

                                 cfg
                             })

        shared_cfg <- plyr::rbind.fill(shared_cfg)

        shared_cfgs <- lapply(share.config, function(dummy) { shared_cfg })
        names(shared_cfgs) <- share.config
        shared_cfgs <- shared_cfgs[setdiff(share.config, names(specific))]

    } else {
        shared_cfg <- NULL
        shared_cfgs <- list()
    }

    ## =============================================================================================
    ## Specific configs
    ## =============================================================================================

    if (length(specific) > 0L) {
        cfgs <- mapply(specific, names(specific),
                       SIMPLIFY = FALSE,
                       FUN = function(config, clus_type)
                       {
                           config_names <- names(config)

                           if (!is.list(config) || is.null(config_names))
                               stop("All parameters must be named lists.")

                           cfg <- mapply(config, config_names, SIMPLIFY = FALSE,
                                         FUN = function(config_args, fun) {
                                             cfg <- do.call(expand.grid,
                                                            enlist(foo = fun,
                                                                   stringsAsFactors = FALSE,
                                                                   dots = config_args))

                                             names(cfg)[1L] <- type

                                             cfg
                                         })

                           cfg <- plyr::rbind.fill(cfg)

                           if (clus_type %in% share.config)
                               cfg <- plyr::rbind.fill(shared_cfg, cfg)

                           ## return
                           cfg
                       })

        cfgs <- c(cfgs, shared_cfgs)

    } else {
        cfgs <- shared_cfgs
    }

    ## return
    cfgs
}

#' Create clustering configurations.
#'
#' Create configurations for [compare_clusterings()]
#'
#' @export
#'
#' @param k A numeric vector with one or more elements specifying the number of clusters to test.
#' @param types Clustering types. It must be any combination of (possibly abbreviated): partitional,
#'   hierarchical, fuzzy, tadpole.
#' @param controls A named list of [tsclust-controls]. `NULL` means defaults. See details.
#' @param preprocs Preprocessing configurations. See details.
#' @param distances Distance configurations. See details.
#' @param centroids Centroid configurations. See details.
#'
#' @details
#'
#' This function uses [base::expand.grid()] and [base::merge()].
#'
#' Preprocessing, distance and centroid configurations are specified with the helper function
#' [pdc_configs()], refer to the examples in [compare_clusterings()] to see how this is used.
#'
#' The controls list may be specified with the usual [tsclust-controls] functions. The names  of
#' the list must correspond to "partitional", "hierarchical", "fuzzy" or "tadpole" clustering.
#' Again, please refer to the examples in [compare_clusterings()].
#'
#' @return
#'
#' A list for each clustering type, each of which includes a data frame with the computed and merged
#' configurations. Each data frame has an extra attribute `num.configs` specifying the number of
#' configurations.
#'
compare_clusterings_configs <- function(types = c("p", "h", "f"), k = 2L, controls = NULL,
                                        preprocs = pdc_configs("preproc", none = list()),
                                        distances = pdc_configs("distance", dtw_basic = list()),
                                        centroids = pdc_configs("centroid", default = list()))
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    types <- match.arg(types, supported_clusterings, TRUE)

    ## ---------------------------------------------------------------------------------------------
    ## Check controls specification
    ## ---------------------------------------------------------------------------------------------

    if (is.null(controls)) {
        controls <- lapply(types, function(type) { do.call(paste0(type, "_control"), list()) })
        names(controls) <- types

    } else if (!is.list(controls) || is.null(names(controls))) {
        stop("The 'controls' argument must be NULL or a named list")

    } else if (!all(types %in% names(controls))) {
        stop("The names of the 'controls' argument do not correspond to the provided 'types'")

    } else {
        controls <- controls[intersect(names(controls), types)]
    }

    ## ---------------------------------------------------------------------------------------------
    ## Check preprocessings specification
    ## ---------------------------------------------------------------------------------------------

    if (missing(preprocs))
        force(preprocs)
    else if (!is.list(preprocs) || is.null(names(preprocs)))
        stop("The 'preprocs' argument must be a named list")
    else if (!all(types %in% names(preprocs)))
        stop("The names of the 'preprocs' argument do not correspond to the provided 'types'")

    preprocs <- preprocs[intersect(names(preprocs), types)]

    ## ---------------------------------------------------------------------------------------------
    ## Check distance specification
    ## ---------------------------------------------------------------------------------------------

    if (missing(distances))
        force(distances)
    else if (!is.list(distances) || is.null(names(distances)))
        stop("The 'distances' argument must be a named list")
    else if (!all(setdiff(types, "tadpole") %in% names(distances)))
        stop("The names of the 'distances' argument do not correspond to the provided 'types'")

    distances <- distances[intersect(names(distances), types)]

    ## ---------------------------------------------------------------------------------------------
    ## Check centroids specification
    ## ---------------------------------------------------------------------------------------------

    if (missing(centroids))
        force(centroids)
    else if (!is.list(centroids) || is.null(names(centroids)))
        stop("The 'centroids' argument must be a named list")
    else if (!all(types %in% names(centroids)))
        stop("The names of the 'centroids' argument do not correspond to the provided 'types'")

    centroids <- centroids[intersect(names(centroids), types)]

    ## =============================================================================================
    ## Create configs
    ## =============================================================================================

    ## return here
    mapply(types,
           controls[types],
           preprocs[types],
           distances[types],
           centroids[types],
           SIMPLIFY = FALSE,
           FUN = function(type, control, preproc, distance, centroid) {
               ## ----------------------------------------------------------------------------------
               ## Check control type
               ## ----------------------------------------------------------------------------------

               expected_control_class <- switch(type,
                                                partitional = "PtCtrl",
                                                hierarchical = "HcCtrl",
                                                fuzzy = "FzCtrl",
                                                tadpole = "TpCtrl")

               if (class(control) != expected_control_class) stop("Invalid ", type, " control")

               ## ----------------------------------------------------------------------------------
               ## Control configs
               ## ----------------------------------------------------------------------------------

               cfg <- switch(type,
                             partitional = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                pam.precompute = control$pam.precompute,
                                                iter.max = control$iter.max,
                                                nrep = control$nrep,
                                                symmetric = control$symmetric,
                                                stringsAsFactors = FALSE))
                             },
                             hierarchical = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                method = list(control$method),
                                                symmetric = control$symmetric,
                                                stringsAsFactors = FALSE))
                             },
                             fuzzy = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                fuzziness = control$fuzziness,
                                                iter.max = control$iter.max,
                                                delta = control$delta,
                                                stringsAsFactors = FALSE))
                             },
                             tadpole = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                dc = list(control$dc),
                                                window.size = control$window.size,
                                                lb = control$lb,
                                                stringsAsFactors = FALSE))
                             })

               ## ----------------------------------------------------------------------------------
               ## Merge configs
               ## ----------------------------------------------------------------------------------

               ## preproc
               if (!is.null(preproc) && nrow(preproc)) {
                   nms <- names(preproc)
                   nms_args <- nms != "preproc"

                   if (any(nms_args))
                       names(preproc)[nms_args] <- paste0(nms[nms_args], "_preproc")

                   cfg <- merge(cfg, preproc, all = TRUE)
               }

               ## distance
               if (type != "tadpole" && !is.null(distance) && nrow(distance)) {
                   nms <- names(distance)
                   nms_args <- nms != "distance"

                   if (any(nms_args))
                       names(distance)[nms_args] <- paste0(nms[nms_args], "_distance")

                   cfg <- merge(cfg, distance, all = TRUE)
               }

               ## centroid
               if (!is.null(centroid) && nrow(centroid)) {
                   nms <- names(centroid)
                   nms_args <- nms != "centroid"

                   if (any(nms_args))
                       names(centroid)[nms_args] <- paste0(nms[nms_args], "_centroid")

                   cfg <- merge(cfg, centroid, all = TRUE)
               }

               ## for info
               attr(cfg, "num.configs") <- switch(type,
                                                  partitional = length(k) * sum(cfg$nrep),
                                                  hierarchical = length(k) * length(cfg$method[[1L]]) * nrow(cfg),
                                                  fuzzy = length(k) * nrow(cfg),
                                                  tadpole = length(k) * length(cfg$dc[[1L]]) * nrow(cfg))
               ## return mapply
               cfg
           })
}

#' Compare different clustering configurations
#'
#' Compare many different clustering algorithms with support for parallelization.
#'
#' @export
#'
#' @param series A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced to a list row-wise.
#' @param types Clustering types. It must be any combination of (possibly abbreviated): partitional,
#'   hierarchical, fuzzy, tadpole.
#' @param ... Further arguments for [tsclust()], `score.clus` or `pick.clus`.
#' @param configs The list of data frames with the desired configurations to run. See
#'   [pdc_configs()] and [compare_clusterings_configs()].
#' @param seed Seed for random reproducibility.
#' @param trace Logical indicating that more output should be printed to screen.
#' @param packages A character vector with the names of any packages needed for any functions used
#'   (distance, centroid, preprocessing, etc.). The name "dtwclust" is added automatically. Relevant
#'   for parallel computation.
#' @param score.clus A function that gets the list of results (and `...`) and scores each one. It
#'   may also be a named list of functions, one for each type of clustering. See Scoring section.
#' @param pick.clus A function that gets the result from `score.clus` as first argument, as well as
#'   the objects returned by [tsclust()] (and elements of `...`) and picks the best result.
#' @param shuffle.configs Randomly shuffle the order of configs, which can be useful to balance load
#'   when using parallel computation.
#' @param return.objects Logical indicating whether the objects from returned by [tsclust()] should
#'   be given in the result.
#'
#' @details
#'
#' This function calls [tsclust()] with different configurations and evaluates the results with the
#' provided functions. Parallel support is included. See the examples.
#'
#' Parameters specified in `configs` whose values are `NA` will be ignored automatically.
#'
#' The scoring and picking functions are for convenience, if they are not specified, the `scores`
#' and `pick` elements of the result will be `NULL`.
#'
#' @return
#'
#' A list with:
#'
#' - `results`: A list of data frames with the flattened configs and the corresponding scores
#'   returned by `score.clus`.
#' - `scores`: The scores given by `score.clus`.
#' - `pick`: The object returned by `pick.clus`.
#' - `proc_time`: The measured execution time, using [base::proc.time()].
#'
#' The cluster objects are also returned if `return.objects` `=` `TRUE`.
#'
#' @section Scoring:
#'
#'   The clustering results are organized in a *list of lists* in the following way (where only
#'   applicable `types` exist; first-level list names in bold):
#'
#'   - **partitional** - list with
#'     + Clustering results from first partitional config
#'     + etc.
#'   - **hierarchical** - list with
#'     + Clustering results from first hierarchical config
#'     + etc.
#'   - **fuzzy** - list with
#'     + Clustering results from first fuzzy config
#'     + etc.
#'   - **tadpole** - list with
#'     + Clustering results from first tadpole config
#'     + etc.
#'
#'   If `score.clus` is a function, it will be applied to the available partitional, hierarchical,
#'   fuzzy and/or tadpole results via:
#'
#'   ```
#'   scores <- lapply(list_of_lists, score.clus, ...)
#'   ```
#'
#'   Otherwise, `score.clus` should be a list of functions with the same names as the list above,
#'   so that `score.clus$partitional` is used to score `list_of_lists$partitional` and so on (via
#'   [base::mapply()] with `SIMPLIFY` `=` `FALSE`).
#'
#'   Therefore, the scores returned shall always be a list of lists with first-level names as above.
#'   These scores and the list of clustering results are given to `pick.clus` as first and second
#'   arguments respectively, followed by `...`.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @seealso
#'
#' [compare_clusterings_configs()], [tsclust()]
#'
#' @example inst/comparison-examples.R
#'
compare_clusterings <- function(series = NULL, types = c("p", "h", "f", "t"), ...,
                                configs = compare_clusterings_configs(types),
                                seed = NULL, trace = FALSE, packages = character(0L),
                                score.clus = function(...) stop("No scoring"),
                                pick.clus = function(...) stop("No picking"),
                                shuffle.configs = FALSE, return.objects = FALSE)
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    tic <- proc.time()
    set.seed(seed)

    if (is.null(series)) stop("No series provided.")
    types <- match.arg(types, supported_clusterings, TRUE)

    ## coerce to list if necessary
    series <- any2list(series)
    check_consistency(series, "vltslist")

    if (!is.function(score.clus) && !(is.list(score.clus) && all(sapply(score.clus, is.function))))
        stop("Invalid evaluation function(s)")
    else if (is.list(score.clus)) {
        if (!all(types %in% names(score.clus)))
            stop("The names of the 'score.clus' argument do not correspond to the provided 'types'")

        score.clus <- score.clus[types]
    }

    if (!is.function(pick.clus)) stop("Invalid pick function")

    ## ---------------------------------------------------------------------------------------------
    ## Misc parameters
    ## ---------------------------------------------------------------------------------------------

    packages <- unique(c("dtwclust", packages))
    dots <- list(...)
    configs <- configs[types]

    if (any(sapply(configs, is.null)))
        stop("The configuration for one of the chosen clustering types is missing.")

    if (shuffle.configs) {
        configs <- lapply(configs, function(config) {
            config[sample(nrow(config)), , drop = FALSE]
        })
    }

    ## ---------------------------------------------------------------------------------------------
    ## Obtain random seeds
    ## ---------------------------------------------------------------------------------------------

    num_seeds <- cumsum(sapply(configs, nrow))
    seeds <- rngtools::RNGseq(num_seeds[length(num_seeds)], seed = seed, simplify = FALSE)
    seeds <- mapply(c(1L, num_seeds[-length(num_seeds)] + 1L),
                    num_seeds,
                    SIMPLIFY = FALSE,
                    FUN = function(first, last) { seeds[first:last] })
    names(seeds) <- names(configs)

    ## =============================================================================================
    ## Preprocessings
    ## =============================================================================================

    if (trace) cat("Preprocessing series...\n")

    processed_series <- lapply(configs, function(config) {
        preproc_cols <- grepl("_?preproc$", names(config))
        preproc_df <- config[, preproc_cols, drop = FALSE]

        if (nrow(preproc_df) > 0L)
            preproc_df <- unique(preproc_df)
        else
            return(NULL)

        preproc_args <- grepl("_preproc$", names(preproc_df))

        lapply(seq_len(nrow(preproc_df)), function(i) {
            if (preproc_df$preproc[i] != "none") {
                preproc_fun <- get_from_callers(preproc_df$preproc[i], "function")

                this_config <- preproc_df[i, , drop = FALSE]
                names(this_config) <- gsub("_preproc", "", names(this_config))

                if (any(preproc_args)) {
                    preproc_args <- as.list(this_config)[preproc_args]
                    preproc_args <- preproc_args[!sapply(preproc_args, is.na)]

                } else
                    preproc_args <- list()

                ret <- do.call(preproc_fun,
                               enlist(series,
                                      dots = preproc_args))

                attr(ret, "config") <- as.list(this_config)

            } else {
                ret <- series
                attr(ret, "config") <- list(preproc = "none")
            }

            ret
        })
    })

    names(processed_series) <- names(configs)

    ## =============================================================================================
    ## Matrix allocation
    ## =============================================================================================

    matrices_allocated <- FALSE
    if (any(types != "tadpole")) {
        allocate_gcm <- any(sapply(setdiff(types, "tadpole"), function(type) {
            "dtw_basic" %in% configs[[type]]$distance && !("gcm" %in% colnames(configs[[type]]))
        }))

        allocate_logs <- any(sapply(setdiff(types, "tadpole"), function(type) {
            "gak" %in% configs[[type]]$distance && !("logs" %in% colnames(configs[[type]]))
        }))

        if (allocate_gcm || allocate_logs)
            N <- max(sapply(processed_series, function(series_by_type) {
                max(sapply(series_by_type, function(series) {
                    max(sapply(series, NROW))
                }))
            }))

        if (allocate_gcm && is.null(dots$gcm)) {
            dots$gcm <- matrix(0, 2L, N + 1L)
            matrices_allocated <- TRUE
        }

        if (allocate_logs && is.null(dots$logs)) {
            dots$logs <- matrix(0, N + 1L, 3L)
            matrices_allocated <- TRUE
        }
    }

    ## =============================================================================================
    ## Clusterings
    ## =============================================================================================

    if (trace) cat("\n")

    objs_by_type <- mapply(configs, names(configs), seeds, SIMPLIFY = FALSE, FUN = function(config, type, seeds) {
        if (trace) message("=================================== Performing ",
                           type,
                           " clusterings ===================================\n")

        series <- processed_series[[type]]

        ## -----------------------------------------------------------------------------------------
        ## distance entries to re-register in parallel workers
        ## -----------------------------------------------------------------------------------------

        dist_names <- unique(config$distance)
        dist_entries <- lapply(dist_names, function(dist) { proxy::pr_DB$get_entry(dist) })
        names(dist_entries) <- dist_names

        ## -----------------------------------------------------------------------------------------
        ## export any necessary preprocessing and centroid functions
        ## -----------------------------------------------------------------------------------------

        custom_preprocs <- setdiff(unique(config$preproc), "none")
        custom_centroids <- setdiff(unique(config$centroid), c("default", centroids_included))

        for (custom_preproc in custom_preprocs)
            assign(custom_preproc, get_from_callers(custom_preproc, "function"))

        for (custom_centroid in custom_centroids)
            assign(custom_centroid, get_from_callers(custom_centroid, "function"))

        export <- c("dots", "trace", "centroids_included", "matrices_allocated",
                    "check_consistency", "enlist", "subset_dots", "get_from_callers",
                    custom_preprocs, custom_centroids)

        ## -----------------------------------------------------------------------------------------
        ## perform clusterings
        ## -----------------------------------------------------------------------------------------

        config_split <- split_parallel(config, 1L)
        seeds_split <- split_parallel(seeds)
        cfg <- NULL # CHECK complains with NOTE regarding undefined global

        objs <- foreach(
            cfg = config_split,
            seed = seeds_split,
            .combine = c,
            .multicombine = TRUE,
            .packages = packages,
            .export = export
        ) %op% {
            chunk <- lapply(seq_len(nrow(cfg)), function(i) {
                if (trace) {
                    message("-------------- Using configuration: --------------")
                    print(cfg[i, , drop = FALSE])
                }

                ## ---------------------------------------------------------------------------------
                ## obtain args from configuration
                ## ---------------------------------------------------------------------------------

                args <- lapply(c("preproc", "distance", "centroid"),
                               function(func) {
                                   col_ids <- grepl(paste0("_", func, "$"), names(cfg))

                                   if (cfg[[func]][i] != "none" && any(col_ids)) {
                                       this_args <- as.list(cfg[i, col_ids, drop = FALSE])
                                       names(this_args) <- sub(paste0("_", func),
                                                               "",
                                                               names(this_args))
                                       ## return
                                       this_args[!sapply(this_args, is.na)]

                                   } else list()
                               })

                names(args) <- c("preproc", "dist", "cent")
                args <- do.call(tsclust_args, args = args)

                ## ---------------------------------------------------------------------------------
                ## controls for this configuration
                ## ---------------------------------------------------------------------------------

                control_fun <- match.fun(paste0(type, "_control"))
                control_args <- subset_dots(as.list(cfg[i, , drop = FALSE]), control_fun)
                control_args <- lapply(control_args, unlist, recursive = FALSE)
                control <- do.call(control_fun, control_args)

                ## ---------------------------------------------------------------------------------
                ## get processed series
                ## ---------------------------------------------------------------------------------

                preproc_char <- cfg$preproc[i]
                this_preproc_config <- c(list(preproc = preproc_char), args$preproc)

                for (this_series in series) {
                    this_series_config <- attr(this_series, "config")
                    if (identical(this_series_config[names(this_preproc_config)], this_preproc_config))
                        break
                }

                ## ---------------------------------------------------------------------------------
                ## distance entry to re-register in parallel worker
                ## ---------------------------------------------------------------------------------

                if (type != "tadpole") {
                    distance <- cfg$distance[i]
                    dist_entry <- dist_entries[[distance]]

                    if (!check_consistency(dist_entry$names[1L], "dist"))
                        do.call(proxy::pr_DB$set_entry, dist_entry)

                } else distance <- "dtw_basic" ## dummy

                ## ---------------------------------------------------------------------------------
                ## centroid for this configuration
                ## ---------------------------------------------------------------------------------

                centroid_char <- cfg$centroid[i]

                ## ---------------------------------------------------------------------------------
                ## call tsclust
                ## ---------------------------------------------------------------------------------

                this_args <- enlist(series = this_series,
                                    type = type,
                                    k = cfg$k[[i]],
                                    distance = distance,
                                    seed = seed[[i]],
                                    trace = trace,
                                    args = args,
                                    control = control,
                                    dots = dots)

                if (type == "tadpole") this_args$distance <- NULL

                if (centroid_char == "default") {
                    ## do not specify centroid
                    tsc <- do.call(tsclust, this_args, quote = TRUE)

                } else if (centroid_char %in% centroids_included) {
                    ## with included centroid
                    tsc <- do.call(tsclust,
                                   enlist(centroid = centroid_char, dots = this_args),
                                   quote = TRUE)

                } else {
                    ## with centroid function
                    delayedAssign("centroid", get_from_callers(centroid_char, "function"))

                    tsc <- do.call(tsclust,
                                   enlist(centroid = centroid, dots = this_args),
                                   quote = TRUE)
                }

                if (inherits(tsc, "TSClusters")) tsc <- list(tsc)

                tsc <- lapply(tsc, function(tsc) {
                    tsc@preproc <- preproc_char

                    if (preproc_char != "none")
                        tsc@family@preproc <- get_from_callers(preproc_char, "function")
                    if (centroid_char != "default")
                        tsc@centroid <- centroid_char

                    if (matrices_allocated) {
                        tsc@dots$gcm <- tsc@dots$logs <- NULL
                        tsc@args <- lapply(tsc@args, function(arg) {
                            arg$gcm <- arg$logs <- NULL
                            arg
                        })
                    }

                    tsc
                })

                tsc
            })

            ## return chunk
            unlist(chunk, recursive = FALSE, use.names = FALSE)
        }

        ## return TSClusters
        objs
    })

    ## =============================================================================================
    ## Evaluations
    ## =============================================================================================

    if (is.function(score.clus))
        scores <- try(lapply(objs_by_type, score.clus, ...), silent = TRUE)
    else
        scores <- try(mapply(objs_by_type, score.clus[names(objs_by_type)],
                             SIMPLIFY = FALSE,
                             MoreArgs = dots,
                             FUN = function(objs, score_fun, ...) { score_fun(objs, ...) }),
                      silent = TRUE)

    if (inherits(scores, "try-error")) {
        warning("The score.clus function(s) did not execute successfully.")
        scores <- NULL
        pick <- NULL

    } else {
        pick <- try(pick.clus(scores, objs_by_type, ...), silent = TRUE)

        if (inherits(pick, "try-error")) {
            warning("The pick.clus function did not execute successfully.")
            pick <- NULL
        }
    }

    ## =============================================================================================
    ## Data frame with results
    ## =============================================================================================

    k <- unlist(configs[[1L]]$k[1L])
    i_cfg <- 1L
    config_ids <- list()

    for (nr in sapply(configs, nrow)) {
        config_ids <- c(config_ids, list(seq(from = i_cfg, by = 1L, length.out = nr)))
        i_cfg <- i_cfg + nr
    }

    configs_out <- mapply(configs, config_ids, types, SIMPLIFY = FALSE, FUN = function(config, ids, type) {
        config <- data.frame(config_id = paste0("config", ids), config)

        dfs <- switch(type,
                      partitional = {
                          if (length(k) > 1L || any(config$nrep > 1L))
                              lapply(seq_len(nrow(config)), function(i) {
                                  this_config <- config[i, , drop = FALSE]
                                  rep <- 1L:this_config$nrep
                                  this_config <- this_config[setdiff(
                                      names(this_config), c("k", "nrep")
                                  )]

                                  df <- expand.grid(rep = rep, k = k)
                                  make_unique_ids(df, this_config) ## see EOF
                              })
                          else
                              config
                      },
                      hierarchical = {
                          if (length(k) > 1L || any(lengths(config$method) > 1L))
                              lapply(seq_len(nrow(config)), function(i) {
                                  this_config <- config[i, , drop = FALSE]
                                  method <- unlist(this_config$method)
                                  this_config <- this_config[setdiff(
                                      names(this_config), c("k", "method")
                                  )]

                                  df <- expand.grid(k = k, method = method)
                                  make_unique_ids(df, this_config) ## see EOF
                              })
                          else
                              config
                      },
                      fuzzy = {
                          if (length(k) > 1L)
                              lapply(seq_len(nrow(config)), function(i) {
                                  this_config <- config[i, , drop = FALSE]
                                  this_config <- this_config[setdiff(
                                      names(this_config), c("k")
                                  )]

                                  df <- expand.grid(k = k)
                                  make_unique_ids(df, this_config) ## see EOF
                              })
                          else
                              config
                      },
                      tadpole = {
                          if (length(k) > 1L || any(lengths(config$dc) > 1L))
                              lapply(seq_len(nrow(config)), function(i) {
                                  this_config <- config[i, , drop = FALSE]
                                  dc <- unlist(this_config$dc)
                                  this_config <- this_config[setdiff(
                                      names(this_config), c("k", "dc")
                                  )]

                                  df <- expand.grid(k = k, dc = dc)
                                  make_unique_ids(df, this_config) ## see EOF
                              })
                          else
                              config
                      })

        ## return mapply
        plyr::rbind.fill(dfs)
    })

    ## ---------------------------------------------------------------------------------------------
    ## Add scores
    ## ---------------------------------------------------------------------------------------------

    ## in case ordering is required below
    if (shuffle.configs)
        configs_cols <- lapply(configs_out, function(config) {
            setdiff(colnames(config), c("config_id", "rep"))
        })

    if (!is.null(scores)) {
        results <- try(mapply(configs_out, scores,
                              SIMPLIFY = FALSE,
                              FUN = function(config, score) {
                                  cbind(config, as.data.frame(score))
                              }),
                       silent = TRUE)

        if (inherits(results, "try-error")) {
            warning("The scores could not be appended to the results data frame.")
            results <- configs_out
        }

    } else {
        results <- configs_out
    }

    ## =============================================================================================
    ## List with all results
    ## =============================================================================================

    results <- list(results = results, scores = scores, pick = pick, proc_time = proc.time() - tic)

    if (return.objects) {
        objs_by_type <- mapply(objs_by_type, results$results,
                               SIMPLIFY = FALSE,
                               FUN = function(objs, res) {
                                   names(objs) <- res$config_id
                                   objs
                               })

        results <- c(results, objects = objs_by_type)
    }

    if (shuffle.configs)
        results$results <- mapply(results$results, configs_cols[names(results$results)],
                                  SIMPLIFY = FALSE,
                                  FUN = function(result, cols) {
                                      order_args <- as.list(result[cols])
                                      names(order_args) <- NULL
                                      result[do.call(order, order_args), , drop = FALSE]
                                  })

    ## return results
    results
}

## =================================================================================================
## compare_clusterings helpers
## =================================================================================================

make_unique_ids <- function(df, this_config) {
    rownames(this_config) <- NULL
    this_config <- cbind(this_config[, 1L, drop = FALSE], df, this_config[, -1L, drop = FALSE])
    nr <- nrow(this_config)
    if (nr > 1L) this_config$config_id <- paste0(this_config$config_id, "_", 1L:nr)
    this_config
}
