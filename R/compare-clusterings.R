#' Helper function for preprocessing/distance/centroid configurations
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

    share_missing <- missing(share.config)
    share.config <- match.arg(share.config, supported_clusterings, TRUE)

    if (type == "distance") {
        if (!is.null(specific$tadpole) || (!share_missing && "tadpole" %in% share.config))
            warning("TADPole ignores distance configurations.")

        specific$tadpole <- NULL
        share.config <- setdiff(share.config, "tadpole")
    }

    ## =============================================================================================
    ## Shared configs
    ## =============================================================================================

    if (length(shared) > 0L && length(share.config) > 0L) {
        shared_names <- names(shared)

        ## careful, singular and plural below
        shared_cfg <- Map(shared, shared_names, f = function(shared_args, fun) {
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
        cfgs <- Map(specific, names(specific), f = function(config, clus_type) {
            config_names <- names(config)

            if (!is.list(config) || is.null(config_names))
                stop("All parameters must be named lists.")

            cfg <- Map(config, config_names, f = function(config_args, fun) {
                cfg <- do.call(expand.grid,
                               enlist(foo = fun,
                                      stringsAsFactors = FALSE,
                                      dots = config_args))

                names(cfg)[1L] <- type
                cfg
            })

            cfg <- plyr::rbind.fill(cfg)
            if (clus_type %in% share.config) cfg <- plyr::rbind.fill(shared_cfg, cfg)
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
    Map(types, controls[types], preprocs[types], distances[types], centroids[types],
        f = function(type, control, preproc, distance, centroid) {
            ## -------------------------------------------------------------------------------------
            ## Control configs
            ## -------------------------------------------------------------------------------------

            if (class(control) != control_classes[type]) stop("Invalid ", type, " control")

            cfg <- switch(
                type,
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
                                   symmetric = control$symmetric,
                                   stringsAsFactors = FALSE))
                },
                tadpole = {
                    do.call(expand.grid,
                            enlist(k = list(k),
                                   dc = list(control$dc),
                                   window.size = control$window.size,
                                   lb = control$lb,
                                   stringsAsFactors = FALSE))
                }
            )

            ## -------------------------------------------------------------------------------------
            ## Merge configs
            ## -------------------------------------------------------------------------------------

            ## preproc
            if (!is.null(preproc) && nrow(preproc)) {
                nms <- names(preproc)
                nms_args <- nms != "preproc"
                if (any(nms_args)) names(preproc)[nms_args] <- paste0(nms[nms_args], "_preproc")
                cfg <- merge(cfg, preproc, all = TRUE)
            }

            ## distance
            if (type != "tadpole" && !is.null(distance) && nrow(distance)) {
                nms <- names(distance)
                nms_args <- nms != "distance"
                if (any(nms_args)) names(distance)[nms_args] <- paste0(nms[nms_args], "_distance")
                cfg <- merge(cfg, distance, all = TRUE)
            }

            ## centroid
            if (!is.null(centroid) && nrow(centroid)) {
                nms <- names(centroid)
                nms_args <- nms != "centroid"
                if (any(nms_args)) names(centroid)[nms_args] <- paste0(nms[nms_args], "_centroid")
                cfg <- merge(cfg, centroid, all = TRUE)
            }

            ## for info
            attr(cfg, "num.configs") <- switch(
                type,
                partitional = length(k) * sum(cfg$nrep),
                hierarchical = length(k) * length(cfg$method[[1L]]) * nrow(cfg),
                fuzzy = length(k) * nrow(cfg),
                tadpole = length(k) * length(cfg$dc[[1L]]) * nrow(cfg)
            )

            ## return Map
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
#' @param score.clus A function that gets the list of results (and `...`) and scores each one. It
#'   may also be a named list of functions, one for each type of clustering. See Scoring section.
#' @param pick.clus A function that to pick the best result. See Picking section.
#' @param shuffle.configs Randomly shuffle the order of configs, which can be useful to balance load
#'   when using parallel computation.
#' @param return.objects Logical indicating whether the objects returned by [tsclust()] should be
#'   given in the result.
#' @param packages A character vector with the names of any packages needed for any functions used
#'   (distance, centroid, preprocessing, etc.). The name "dtwclust" is added automatically. Relevant
#'   for parallel computation.
#' @param .errorhandling This will be passed to [foreach::foreach()]. See Parallel section below.
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
#' - `seeds`: A list of lists with the random seeds computed for each configuration.
#'
#' The cluster objects are also returned if `return.objects` `=` `TRUE`.
#'
#' @section Parallel computation:
#'
#'   The configurations for each clustering type can be evaluated in parallel with the \pkg{foreach}
#'   package. A parallel backend can be registered, e.g., with \pkg{doParallel}.
#'
#'   If the `.errorhandling` parameter is changed to "pass" and a custom `score.clus` function is
#'   used, said function should be able to deal with possible error objects.
#'
#'   If it is changed to "remove", it might not be possible to attach the scores to the results data
#'   frame, or it may be inconsistent. Additionally, if `return.objects` is `TRUE`, the names given
#'   to the objects might also be inconsistent.
#'
#'   Parallelization can incur a lot of deep copies of data when returning the cluster objects,
#'   since each one will likely contain a copy of `datalist`. If you want to avoid this, consider
#'   specifying `score.clus` and setting `return.objects` to `FALSE`.
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
#'   [base::Map()]).
#'
#'   Therefore, the scores returned shall always be a list of lists with first-level names as above.
#'
#' @section Picking:
#'
#'   If `return.objects` is `TRUE`, the scores and the list of clustering results are given to
#'   `pick.clus` as first and second arguments respectively, followed by `...`. Otherwise,
#'   `pick.clus` will receive only the scores and the contents of `...` (since the results will not
#'   be returned by the preceding step).
#'
#' @section Limitations:
#'
#'   Note that the configurations returned by the helper functions assign special names to
#'   preprocessing/distance/centroid arguments, and these names are used internally to recognize
#'   them.
#'
#'   If some of these arguments are more complicated (e.g. matrices) and should *not* be expanded,
#'   consider passing them directly via the ellipsis (`...`) instead of using [pdc_configs()]. This
#'   assumes that said arguments can be passed to all functions without affecting their results.
#'
#'   The distance matrices (if calculated) are not re-used across configurations. Given the way the
#'   configurations are created, this shouldn't matter, because clusterings with arguments that can
#'   use the same distance matrix are already grouped together by [compare_clusterings_configs()]
#'   and [pdc_configs()].
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
                                seed = NULL, trace = FALSE,
                                score.clus = function(...) stop("No scoring"),
                                pick.clus = function(...) stop("No picking"),
                                shuffle.configs = FALSE, return.objects = FALSE,
                                packages = character(0L), .errorhandling = "stop")
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    tic <- proc.time()
    set.seed(seed)

    score_missing <- missing(score.clus)
    pick_missing <- missing(pick.clus)

    if (is.null(series)) stop("No series provided.")
    types <- match.arg(types, supported_clusterings, TRUE)
    .errorhandling <- match.arg(.errorhandling, c("stop", "remove", "pass"))
    if (!return.objects && score_missing)
        stop("Returning no objects and specifying no scoring function would return no useful results.")

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
    seeds <- Map(c(1L, num_seeds[-length(num_seeds)] + 1L), num_seeds,
                 f = function(first, last) { seeds[first:last] })
    setnames_inplace(seeds, names(configs))

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
                names(this_config) <- sub("_preproc$", "", names(this_config))

                if (any(preproc_args)) {
                    preproc_args <- as.list(this_config)[preproc_args]
                    preproc_args <- preproc_args[!sapply(preproc_args, is.na)]

                } else
                    preproc_args <- list()

                ret <- do.call(preproc_fun,
                               enlist(series,
                                      dots = preproc_args))

                attr(ret, "config") <- as.list(this_config) ## leave version with possible NAs here!

            } else {
                ret <- series
            }

            ret
        })
    })

    ## utils.R
    setnames_inplace(processed_series, names(configs))

    ## =============================================================================================
    ## Matrix allocation
    ## =============================================================================================

    allocate_gcm <- allocate_logs <- allocate_dba <- FALSE
    if (any(types != "tadpole")) {
        allocate_gcm <- any(sapply(setdiff(types, "tadpole"), function(type) {
            any(grepl("^dtw_basic$|^dtw_lb$", configs[[type]]$distance, ignore.case = TRUE)) &&
                !("gcm_distance" %in% colnames(configs[[type]]))
        }))

        allocate_logs <- any(sapply(setdiff(types, "tadpole"), function(type) {
            any(grepl("^gak$", configs[[type]]$distance, ignore.case = TRUE)) &&
                !("logs_distance" %in% colnames(configs[[type]]))
        }))

        allocate_dba <- any(sapply(types, function(type) {
            any(grepl("^dba$", configs[[type]]$centroid, ignore.case = TRUE)) &&
                !("gcm_centroid" %in% colnames(configs[[type]]))
        }))

        if (allocate_gcm || allocate_logs || allocate_dba)
            N <- max(sapply(processed_series, function(series_by_type) {
                max(sapply(series_by_type, function(series) {
                    max(sapply(series, NROW))
                }))
            }))

        if (is.null(dots$gcm)) {
            if (allocate_dba)
                dots$gcm <- matrix(0, N + 1L, N + 1L) ## will also work if dtw_basic + DBA is used
            else if (allocate_gcm)
                dots$gcm <- matrix(0, 2L, N + 1L)

        } else allocate_dba <- allocate_gcm <- FALSE

        if (allocate_logs && is.null(dots$logs))
            dots$logs <- matrix(0, N + 1L, 3L)
        else
            allocate_logs <- FALSE
    }

    ## =============================================================================================
    ## Clusterings
    ## =============================================================================================

    if (trace) cat("\n")

    objs_by_type <- Map(configs, names(configs), seeds, f = function(config, type, seeds) {
        if (trace) message("=================================== Performing ",
                           type,
                           " clusterings ===================================\n")

        series <- processed_series[[type]]

        ## -----------------------------------------------------------------------------------------
        ## distance entries to re-register in parallel workers
        ## -----------------------------------------------------------------------------------------

        if (type != "tadpole") {
            dist_names <- unique(config$distance)
            dist_entries <- lapply(dist_names, function(dist) { proxy::pr_DB$get_entry(dist) })
            setnames_inplace(dist_entries, dist_names)
        }

        ## -----------------------------------------------------------------------------------------
        ## export any necessary preprocessing and centroid functions
        ## -----------------------------------------------------------------------------------------

        custom_preprocs <- setdiff(unique(config$preproc), "none")
        custom_centroids <- setdiff(unique(config$centroid), c("default", centroids_included))

        for (custom_preproc in custom_preprocs)
            assign(custom_preproc, get_from_callers(custom_preproc, "function"))

        for (custom_centroid in custom_centroids)
            assign(custom_centroid, get_from_callers(custom_centroid, "function"))

        export <- c("trace", "score.clus", "return.objects",
                    "dots", "allocate_gcm", "allocate_logs", "allocate_dba",
                    "centroids_included",
                    "check_consistency", "enlist", "subset_dots", "get_from_callers", "setnames_inplace",
                    custom_preprocs, custom_centroids)

        ## -----------------------------------------------------------------------------------------
        ## perform clusterings
        ## -----------------------------------------------------------------------------------------

        force(seeds)
        i <- nrow(config)

        objs <- foreach::foreach(
            i = seq_len(i),
            .combine = list,
            .multicombine = TRUE,
            .packages = packages,
            .export = export,
            .errorhandling = .errorhandling,
            .maxcombine = if (isTRUE(i > 100L)) i else 100L
        ) %op% {
            cfg <- config[i, , drop = FALSE]
            seed <- seeds[[i]]

            if (trace) {
                message("-------------- Using configuration: --------------")
                print(cfg)
                cat("\n")
            }

            ## ---------------------------------------------------------------------------------
            ## obtain args from configuration
            ## ---------------------------------------------------------------------------------

            args <- lapply(c("preproc", "distance", "centroid"),
                           function(func) {
                               col_ids <- grepl(paste0("_", func, "$"), names(cfg))

                               if (cfg[[func]] != "none" && any(col_ids)) {
                                   this_args <- as.list(cfg[, col_ids, drop = FALSE])
                                   names(this_args) <- sub(paste0("_", func, "$"),
                                                           "",
                                                           names(this_args))
                                   ## return
                                   this_args[!sapply(this_args, is.na)]

                               } else list()
                           })

            setnames_inplace(args, c("preproc", "dist", "cent"))
            args <- do.call(tsclust_args, args = args)

            ## ---------------------------------------------------------------------------------
            ## controls for this configuration
            ## ---------------------------------------------------------------------------------

            control_fun <- match.fun(paste0(type, "_control"))
            control_args <- subset_dots(as.list(cfg), control_fun)
            control_args <- lapply(control_args, unlist, recursive = FALSE)
            control <- do.call(control_fun, control_args)

            ## ---------------------------------------------------------------------------------
            ## get processed series
            ## ---------------------------------------------------------------------------------

            preproc_char <- cfg$preproc
            this_preproc_config <- c(list(preproc = preproc_char), args$preproc)

            for (this_series in series) {
                this_series_config <- attr(this_series, "config")
                if (preproc_char == "none" && is.null(this_series_config))
                    break
                if (identical(this_series_config[names(this_preproc_config)], this_preproc_config))
                    break
            }

            ## ---------------------------------------------------------------------------------
            ## distance entry to re-register in parallel worker
            ## ---------------------------------------------------------------------------------

            if (type != "tadpole") {
                distance <- cfg$distance
                dist_entry <- dist_entries[[distance]]

                if (!check_consistency(dist_entry$names[1L], "dist"))
                    do.call(proxy::pr_DB$set_entry, dist_entry)

            } else distance <- NULL ## dummy

            ## ---------------------------------------------------------------------------------
            ## centroid for this configuration
            ## ---------------------------------------------------------------------------------

            centroid_char <- cfg$centroid

            ## ---------------------------------------------------------------------------------
            ## call tsclust
            ## ---------------------------------------------------------------------------------

            this_args <- enlist(series = this_series,
                                type = type,
                                k = unlist(cfg$k),
                                distance = distance,
                                seed = seed,
                                trace = trace,
                                args = args,
                                control = control,
                                error.check = FALSE,
                                dots = dots)

            if (type == "tadpole") this_args$distance <- NULL

            if (centroid_char == "default") {
                ## do not specify centroid
                tsc <- do.call(tsclust, this_args)

            } else if (centroid_char %in% centroids_included) {
                ## with included centroid
                tsc <- do.call(tsclust,
                               enlist(centroid = centroid_char, dots = this_args))

            } else {
                ## with centroid function
                tsc <- do.call(tsclust,
                               enlist(centroid = get_from_callers(centroid_char, "function"),
                                      dots = this_args))
            }

            if (inherits(tsc, "TSClusters")) tsc <- list(tsc)

            ret <- lapply(tsc, function(tsc) {
                tsc@preproc <- preproc_char

                if (preproc_char != "none")
                    tsc@family@preproc <- get_from_callers(preproc_char, "function")
                if (centroid_char != "default")
                    tsc@centroid <- centroid_char

                if (allocate_gcm || allocate_dba) {
                    tsc@dots$gcm <- NULL
                    tsc@args <- lapply(tsc@args, function(arg) {
                        arg$gcm <- NULL
                        arg
                    })
                }

                if (allocate_logs) {
                    tsc@dots$logs <- NULL
                    tsc@args <- lapply(tsc@args, function(arg) {
                        arg$logs <- NULL
                        arg
                    })
                }

                tsc
            })

            ## ---------------------------------------------------------------------------------
            ## evaluate
            ## ---------------------------------------------------------------------------------

            if (!return.objects) {
                if (!is.function(score.clus)) score.clus <- score.clus[[type]]
                ret <- list(do.call(score.clus, enlist(ret, dots = dots)))
            }

            ## return config result from foreach()
            ret
        }

        ## foreach() with one 'iteration' does NOT call .combine function
        if (i == 1L) objs <- list(objs)

        if (.errorhandling == "pass") {
            failed_cfgs <- sapply(objs, function(obj) { inherits(obj, "error") })
            if (any(failed_cfgs)) {
                warning("At least one of the ", type, " configurations resulted in an error.")
                ## make one more list level so unlist() below leaves errors as one object
                objs[failed_cfgs] <- lapply(objs[failed_cfgs], list)
            }
        }

        ## return scores/TSClusters from Map()
        unlist(objs, recursive = FALSE)
    })

    ## =============================================================================================
    ## Evaluations
    ## =============================================================================================

    if (return.objects) {
        if (is.function(score.clus))
            scores <- try(lapply(objs_by_type, score.clus, ...), silent = TRUE)
        else
            scores <- try(mapply(objs_by_type, score.clus[names(objs_by_type)],
                                 SIMPLIFY = FALSE,
                                 MoreArgs = dots,
                                 FUN = function(objs, score_fun, ...) { score_fun(objs, ...) }),
                          silent = TRUE)

        if (inherits(scores, "try-error")) {
            if (!score_missing) warning("The score.clus function(s) did not execute successfully.")
            scores <- NULL
            pick <- NULL

        } else {
            pick <- try(pick.clus(scores, objs_by_type, ...), silent = TRUE)

            if (inherits(pick, "try-error")) {
                if (!pick_missing) warning("The pick.clus function did not execute successfully.")
                pick <- NULL
            }
        }
    } else {
        scores <- lapply(objs_by_type, function(objs) {
            failed_cfgs <- sapply(objs, function(obj) { inherits(obj, "error") })

            if (any(failed_cfgs))
                passed_objs <- objs[!failed_cfgs]
            else
                passed_objs <- objs

            if (length(passed_objs) == 0L) return(objs)

            if (any(sapply(passed_objs, function(score) { is.null(dim(score)) })))
                unlist(passed_objs, recursive = FALSE)
            else
                plyr::rbind.fill(lapply(passed_objs, base::as.data.frame))
        })

        pick <- try(pick.clus(scores, ...), silent = TRUE)
        if (inherits(pick, "try-error")) {
            if (!pick_missing) warning("The pick.clus function did not execute successfully.")
            pick <- NULL
        }
    }

    ## =============================================================================================
    ## Data frame with results
    ## =============================================================================================

    ## create initial IDs
    i_cfg <- 1L
    config_ids <- lapply(sapply(configs, nrow), function(nr) {
        ids <- seq(from = i_cfg, by = 1L, length.out = nr)
        i_cfg <<- i_cfg + nr
        ids
    })

    ## change to config names and assign to seeds (before flattening)
    config_ids <- Map(config_ids, seeds, f = function(ids, seed) {
        nms <- paste0("config", ids)
        try(setnames_inplace(seed, nms), silent = TRUE)
        nms
    })

    ## flatten
    configs_out <- Map(configs, config_ids, types, f = function(config, ids, type) {
        config <- data.frame(config_id = ids, config, stringsAsFactors = FALSE)
        k <- unlist(config$k[1L])

        dfs <- switch(
            type,
            partitional = {
                if (length(k) > 1L || any(config$nrep > 1L))
                    lapply(seq_len(nrow(config)), function(i) {
                        this_config <- config[i, , drop = FALSE]
                        rep <- 1L:this_config$nrep
                        this_config <- this_config[setdiff(names(this_config), c("k", "nrep"))]
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
                        this_config <- this_config[setdiff(names(this_config), c("k", "method"))]
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
                        this_config <- this_config[setdiff(names(this_config), c("k"))]
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
                        this_config <- this_config[setdiff(names(this_config), c("k", "dc"))]
                        df <- expand.grid(k = k, dc = dc)
                        make_unique_ids(df, this_config) ## see EOF
                    })
                else
                    config
            }
        )

        ## return Map
        plyr::rbind.fill(dfs)
    })

    ## ---------------------------------------------------------------------------------------------
    ## Add scores
    ## ---------------------------------------------------------------------------------------------

    ## in case ordering is required below
    if (shuffle.configs)
        configs_cols <- lapply(configs_out, function(config) {
            setdiff(colnames(config), c("config_id", "rep", "k", "method", "dc",  "window.size", "lb"))
        })

    if (!is.null(scores)) {
        results <- try(Map(configs_out, scores,
                           f = function(config, score) {
                               cbind(config, base::as.data.frame(score))
                           }),
                       silent = TRUE)

        if (inherits(results, "try-error")) {
            warning("The scores could not be appended to the results data frame.")
            results <- configs_out
        }

    } else results <- configs_out

    ## =============================================================================================
    ## List with all results
    ## =============================================================================================

    results <- list(results = results,
                    scores = scores,
                    pick = pick,
                    proc_time = proc.time() - tic,
                    seeds = seeds)

    if (return.objects) {
        setnames_res <- try(Map(objs_by_type, results$results,
                                f = function(objs, res) { setnames_inplace(objs, res$config_id) }),
                            silent = TRUE)

        if (inherits(setnames_res, "try-error"))
            warning("Could not assign names to returned objects.")

        results <- c(results, objects = objs_by_type)
    }

    if (shuffle.configs)
        results$results <- Map(results$results, configs_cols[names(results$results)],
                               f = function(result, cols) {
                                   order_args <- as.list(result[cols])
                                   names(order_args) <- NULL
                                   result[do.call(base::order, order_args), , drop = FALSE]
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
