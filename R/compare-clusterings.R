#' Helper functions for \code{\link{compare_clusterings_configs}}
#'
#' Create preprocessing, distance and centroid configurations for
#' \code{\link{compare_clusterings_configs}}. All functions are based on
#' \code{\link[base]{expand.grid}}.
#'
#' @name config-helpers
#' @rdname config-helpers
#'
#' @param ... Any number of named nested lists with functions and arguments that will be shared by
#'   all clusterings. See details.
#' @param partitional A named list of lists with functions and arguments for partitional
#'   clusterings.
#' @param hierarchical A named list of lists with functions and arguments for hierarchical
#'   clusterings.
#' @param fuzzy A named list of lists with functions and arguments for fuzzy clusterings.
#' @param tadpole A named list of lists with functions and arguments for TADPole clusterings.
#' @param share.config A character vector specifying which clusterings should include the shared
#'   lists (the ones specified in \code{...}). It must be one of (possibly abbreviated):
#'   partitional, hierarchical, fuzzy, tadpole.
#'
#' @details
#'
#' The named lists are interpreted in the following way: the name of each element of the list will
#' be considered to be a function name, and the elements of the nested list will be the possible
#' parameters for the function. Each function must have at least an empty list. The parameters may
#' be vectors that specify different values to be tested.
#'
#' For preprocessing, the special name \code{none} signifies no preprocessing.
#'
#' For centroids, the special name \code{default} leaves the centroid unspecified.
#'
#' Please see the example in \code{\link{compare_clusterings}} to see how this is used.
#'
#' @return
#'
#' A list for each clustering type, each of which includes a data frame with the computed
#' configurations.
#'
NULL

#' @rdname config-helpers
#'
#' @export
#'
preproc_configs <- function(..., partitional, hierarchical, fuzzy, tadpole,
                            share.config = c("p", "h", "f", "t"))
{
    share.config <- match.arg(share.config,
                              c("partitional", "hierarchical", "fuzzy", "tadpole"),
                              TRUE)

    this_call <- match.call(expand.dots = FALSE)
    this_call <- this_call[nzchar(names(this_call))]
    this_call$share.config <- NULL
    this_call$`...` <- NULL

    ## =============================================================================================
    ## Shared configs
    ## =============================================================================================

    shared <- list(...)

    if (length(shared) > 0L) {
        shared_names <- names(shared)

        shared_cfgs <- mapply(shared, shared_names, SIMPLIFY = FALSE,
                              FUN = function(shared_args, preproc) {
                                  do.call(expand.grid,
                                          enlist(preproc = preproc,
                                                 stringsAsFactors = FALSE,
                                                 dots = shared_args))
                              })

        shared_cfgs <- plyr::rbind.fill(shared_cfgs)

        shared_cfgs <- list(partitional = shared_cfgs,
                            hierarchical = shared_cfgs,
                            fuzzy = shared_cfgs,
                            tadpole = shared_cfgs)

    } else {
        shared_cfgs <- list(partitional = NULL,
                            hierarchical = NULL,
                            fuzzy = NULL,
                            tadpole = NULL)
    }

    ## =============================================================================================
    ## Specific configs
    ## =============================================================================================

    if (length(this_call) > 0L) {
        preproc_cfgs <- mapply(this_call, names(this_call),
                               SIMPLIFY = FALSE,
                               FUN = function(config, type)
                               {
                                   config <- eval(config)
                                   config_names <- names(config)

                                   if (!is.list(config) || is.null(config_names))
                                       stop("All parameters must be named lists.")

                                   cfg <- mapply(config, config_names, SIMPLIFY = FALSE,
                                                 FUN = function(config_args, preproc) {
                                                     do.call(expand.grid,
                                                             enlist(preproc = preproc,
                                                                    stringsAsFactors = FALSE,
                                                                    dots = config_args))
                                                 })

                                   cfg <- plyr::rbind.fill(cfg)

                                   if (type %in% share.config)
                                       cfg <- plyr::rbind.fill(shared_cfgs[[type]], cfg)

                                   ## return
                                   cfg
                               })

        names(preproc_cfgs) <- names(this_call)

        shared_cfgs[names(preproc_cfgs)] <- preproc_cfgs
    }

    ## return
    shared_cfgs
}

#' @rdname config-helpers
#'
#' @export
#'
distance_configs <- function(..., partitional, hierarchical, fuzzy,
                             share.config = c("p", "h", "f"))
{
    share.config <- match.arg(share.config,
                              c("partitional", "hierarchical", "fuzzy"),
                              TRUE)

    this_call <- match.call(expand.dots = FALSE)
    this_call <- this_call[nzchar(names(this_call))]
    this_call$share.config <- NULL
    this_call$`...` <- NULL

    ## =============================================================================================
    ## Shared configs
    ## =============================================================================================

    shared <- list(...)

    if (length(shared) > 0L) {
        shared_names <- names(shared)

        shared_cfgs <- mapply(shared, shared_names, SIMPLIFY = FALSE,
                              FUN = function(shared_args, distance) {
                                  do.call(expand.grid,
                                          enlist(distance = distance,
                                                 stringsAsFactors = FALSE,
                                                 dots = shared_args))
                              })

        shared_cfgs <- plyr::rbind.fill(shared_cfgs)

        shared_cfgs <- list(partitional = shared_cfgs,
                            hierarchical = shared_cfgs,
                            fuzzy = shared_cfgs)

    } else {
        shared_cfgs <- list(partitional = NULL,
                            hierarchical = NULL,
                            fuzzy = NULL)
    }

    ## =============================================================================================
    ## Specific configs
    ## =============================================================================================

    if (length(this_call) > 0L) {
        dist_cfgs <- mapply(this_call, names(this_call),
                            SIMPLIFY = FALSE,
                            FUN = function(config, type)
                            {
                                config <- eval(config)
                                config_names <- names(config)

                                if (!is.list(config) || is.null(config_names))
                                    stop("All parameters must be named lists.")

                                cfg <- mapply(config, config_names, SIMPLIFY = FALSE,
                                              FUN = function(config_args, distance) {
                                                  do.call(expand.grid,
                                                          enlist(distance = distance,
                                                                 stringsAsFactors = FALSE,
                                                                 dots = config_args))
                                              })

                                cfg <- plyr::rbind.fill(cfg)

                                if (type %in% share.config)
                                    cfg <- plyr::rbind.fill(shared_cfgs[[type]], cfg)

                                ## return
                                cfg
                            })

        names(dist_cfgs) <- names(this_call)

        shared_cfgs[names(dist_cfgs)] <- dist_cfgs
    }

    ## return
    shared_cfgs
}

#' @rdname config-helpers
#'
#' @export
#'
centroid_configs <- function(..., partitional, hierarchical, fuzzy, tadpole,
                             share.config = c("p", "h", "f", "t"))
{
    share.config <- match.arg(share.config,
                              c("partitional", "hierarchical", "fuzzy", "tadpole"),
                              TRUE)

    this_call <- match.call(expand.dots = FALSE)
    this_call <- this_call[nzchar(names(this_call))]
    this_call$share.config <- NULL
    this_call$`...` <- NULL

    ## =============================================================================================
    ## Shared configs
    ## =============================================================================================

    shared <- list(...)

    if (length(shared) > 0L) {
        shared_names <- names(shared)

        shared_cfgs <- mapply(shared, shared_names, SIMPLIFY = FALSE,
                              FUN = function(shared_args, centroid) {
                                  do.call(expand.grid,
                                          enlist(centroid = centroid,
                                                 stringsAsFactors = FALSE,
                                                 dots = shared_args))
                              })

        shared_cfgs <- plyr::rbind.fill(shared_cfgs)

        shared_cfgs <- list(partitional = shared_cfgs,
                            hierarchical = shared_cfgs,
                            fuzzy = shared_cfgs,
                            tadpole = shared_cfgs)

    } else {
        shared_cfgs <- list(partitional = NULL,
                            hierarchical = NULL,
                            fuzzy = NULL,
                            tadpole = NULL)
    }

    ## =============================================================================================
    ## Specific configs
    ## =============================================================================================

    if (length(this_call) > 0L) {
        centroid_cfgs <- mapply(this_call, names(this_call),
                                SIMPLIFY = FALSE,
                                FUN = function(config, type)
                                {
                                    config <- eval(config)
                                    config_names <- names(config)

                                    if (!is.list(config) || is.null(config_names))
                                        stop("All parameters must be named lists.")

                                    cfg <- mapply(config, config_names, SIMPLIFY = FALSE,
                                                  FUN = function(config_args, centroid) {
                                                      do.call(expand.grid,
                                                              enlist(centroid = centroid,
                                                                     stringsAsFactors = FALSE,
                                                                     dots = config_args))
                                                  })

                                    cfg <- plyr::rbind.fill(cfg)

                                    if (type %in% share.config)
                                        cfg <- plyr::rbind.fill(shared_cfgs[[type]], cfg)

                                    ## return
                                    cfg
                                })

        names(centroid_cfgs) <- names(this_call)

        shared_cfgs[names(centroid_cfgs)] <- centroid_cfgs

    }

    ## return
    shared_cfgs
}

#' Create configurations for \code{\link{compare_clusterings}}
#'
#' Create clustering configurations.
#'
#' @export
#'
#' @param k A numeric vector with one or more elements specifying the number of clusters to test.
#' @param types Clustering types. It must be one of (possibly abbreviated): partitional,
#'   hierarchical, fuzzy, tadpole.
#' @param controls A list of \code{\link{tsclust-controls}}. \code{NULL} means defaults. See
#'   details.
#' @param preprocs Preprocessing configurations. See details.
#' @param distances Distance configurations. See details.
#' @param centroids Centroid configurations. See details.
#'
#' @details
#'
#' This function is based on \code{\link[base]{expand.grid}} and \code{\link[base]{merge}}.
#'
#' Preprocessing, distance and centroid configurations are specified with the helper functions in
#' \code{\link{config-helpers}}, but using the example in \code{\link{compare_clusterings}} as basis
#' might be easier to understand.
#'
#' The controls list must have one nested list of controls for the desired clustering types. Each
#' nested list may be specified with the usual \code{\link{tsclust-controls}} functions. Again,
#' please refer to the example in \code{\link{compare_clusterings}}.
#'
#' @return
#'
#' A list for each clustering type, each of which includes a data frame with the computed and
#' merged configurations. Each data frame has an extra attribute \code{num.configs} specifying the
#' number of configurations.
#'
compare_clusterings_configs <- function(types = c("p", "h", "f", "t"), k = 2L, controls = NULL,
                                        preprocs = preproc_configs(),
                                        distances = distance_configs(),
                                        centroids = centroid_configs())
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    types <- match.arg(types, c("partitional", "hierarchical", "fuzzy", "tadpole"), TRUE)

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
    }

    ## ---------------------------------------------------------------------------------------------
    ## Check preprocessings specification
    ## ---------------------------------------------------------------------------------------------

    if (missing(preprocs)) {
        force(preprocs)
        preprocs <- preprocs[intersect(names(preprocs), types)]

    } else if (!is.list(preprocs) || is.null(names(preprocs))) {
        stop("The 'preprocs' argument must be NULL or a named list")

    } else if (!all(types %in% names(preprocs))) {
        stop("The names of the 'preprocs' argument do not correspond to the provided 'types'")
    }

    ## ---------------------------------------------------------------------------------------------
    ## Check distance specification
    ## ---------------------------------------------------------------------------------------------

    if (missing(distances)) {
        force(distances)
        distances <- distances[intersect(names(distances), types)]

    } else if (!is.list(distances) || is.null(names(distances))) {
        stop("The 'distances' argument must be NULL or a named list")

    } else if (!all(setdiff(types, "tadpole") %in% names(distances))) {
        stop("The names of the 'distances' argument do not correspond to the provided 'types'")
    }

    ## ---------------------------------------------------------------------------------------------
    ## Check centroids specification
    ## ---------------------------------------------------------------------------------------------

    if (missing(centroids)) {
        force(centroids)
        centroids <- centroids[intersect(names(centroids), types)]

    } else if (!is.list(centroids) || is.null(names(centroids))) {
        stop("The 'centroids' argument must be NULL or a named list")

    } else if (!all(types %in% names(centroids))) {
        stop("The names of the 'centroids' argument do not correspond to the provided 'types'")
    }

    ## =============================================================================================
    ## Create configs
    ## =============================================================================================

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

               if (class(control) != expected_control_class)
                   stop("Invalid ", type, " control")

               ## ----------------------------------------------------------------------------------
               ## Control configs
               ## ----------------------------------------------------------------------------------

               cfg <- switch(type,

                             ## --------------------------------------------------------------------
                             ## Partitional configs
                             ## --------------------------------------------------------------------

                             partitional = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                pam.precompute = control$pam.precompute,
                                                iter.max = control$iter.max,
                                                nrep = control$nrep,
                                                stringsAsFactors = FALSE))
                             },

                             ## --------------------------------------------------------------------
                             ## Hierarchical configs
                             ## --------------------------------------------------------------------

                             hierarchical = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                method = list(control$method),
                                                stringsAsFactors = FALSE))
                             },

                             ## --------------------------------------------------------------------
                             ## Fuzzy configs
                             ## --------------------------------------------------------------------

                             fuzzy = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                fuzziness = control$fuzziness,
                                                iter.max = control$iter.max,
                                                delta = control$delta,
                                                stringsAsFactors = FALSE))
                             },

                             ## --------------------------------------------------------------------
                             ## TADPole configs
                             ## --------------------------------------------------------------------

                             tadpole = {
                                 do.call(expand.grid,
                                         enlist(k = list(k),
                                                dc = control$dc,
                                                window.size = control$window.size,
                                                lb = control$lb,
                                                stringsAsFactors = FALSE))
                             })

               ## ----------------------------------------------------------------------------------
               ## Merge configs
               ## ----------------------------------------------------------------------------------

               if (!is.null(preproc) && nrow(preproc)) {
                   nms <- names(preproc)
                   nms_args <- nms != "preproc"

                   if (any(nms_args))
                       names(preproc)[nms_args] <- paste0(nms[nms_args], "_preproc")

                   cfg <- merge(cfg, preproc, all = TRUE)
               }

               if (type != "tadpole" && !is.null(distance) && nrow(distance)) {
                   nms <- names(distance)
                   nms_args <- nms != "distance"

                   if (any(nms_args))
                       names(distance)[nms_args] <- paste0(nms[nms_args], "_distance")

                   cfg <- merge(cfg, distance, all = TRUE)
               }

               if (!is.null(centroid) && nrow(centroid)) {
                   nms <- names(centroid)
                   nms_args <- nms != "centroid"

                   if (any(nms_args))
                       names(centroid)[nms_args] <- paste0(nms[nms_args], "_centroid")

                   cfg <- merge(cfg, centroid, all = TRUE)
               }

               attr(cfg, "num.configs") <- switch(type,
                                                  partitional = length(k) * sum(cfg$nrep),
                                                  hierarchical = length(k) * length(cfg$method[[1L]]) * nrow(cfg),
                                                  fuzzy =, tadpole = length(k) * nrow(cfg))
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
#' @param types Clustering types. It must be one of (possibly abbreviated): partitional,
#'   hierarchical, fuzzy, tadpole.
#' @param ... Further arguments for \code{\link{tsclust}}, \code{score.clus} or \code{pick.clus}.
#' @param configs The list of data frames with the desired configurations to run. See
#'   \code{\link{compare_clusterings_configs}} and \code{\link{config-helpers}}.
#' @param seed Seed for random reproducibility.
#' @param trace Logical indicating that more output should be printed to screen.
#' @param packages A character vector with the names of any packages needed for any functions used
#'   (distance, centroid, preprocessing, etc.). The name "dtwclust" is added automatically. Relevant
#'   for parallel computation.
#' @param score.clus A function that gets the list of results (and \code{...}) and scores each one.
#'   It may also be a named list of functions, one for each type of clustering.
#' @param pick.clus A function that gets the result from \code{score.clus} as first argument, as
#'   well as the objects returned by \code{\link{tsclust}} (and elements of \code{...}) and picks
#'   the best result.
#' @param return.objects Logical indicating whether the objects from returned by
#'   \code{\link{tsclust}} should be given in the result.
#'
#' @details
#'
#' This function calls \code{\link{tsclust}} with different configurations and evaluates the results
#' with the provided functions. Parallel support is included. See the examples.
#'
#' Parameters specified in \code{configs} whose values are \code{NA} will be ignored automatically.
#'
#' The scoring and picking functions are for convenience, if they are not specified, the
#' \code{scores} and \code{pick} elements of the result will be \code{NULL}.
#'
#' @return
#'
#' A list with: \itemize{
#'   \item \code{scores}: The scores given by \code{score.clus}.
#'   \item \code{pick}: The object returned by \code{pick.clus}.
#'   \item \code{proc_time}: The measured function time, using \code{\link[base]{proc.time}}.
#' }
#'
#' The cluster objects are also returned if \code{return.objects} \code{=} \code{TRUE}.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @seealso
#'
#' \code{\link{compare_clusterings_configs}}, \code{\link{tsclust}}
#'
#' @example inst/comparison-examples.R
#'
compare_clusterings <- function(series = NULL, types = c("p", "h", "f", "t"), ...,
                                configs = compare_clusterings_configs(types),
                                seed = NULL, trace = FALSE, packages = character(0L),
                                score.clus = function(...) { stop("No scoring") },
                                pick.clus = function(...) { stop("No picking") },
                                return.objects = FALSE)
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    set.seed(seed)

    tic <- proc.time()

    if (is.null(series)) stop("No data provided")

    types <- match.arg(types, c("partitional", "hierarchical", "tadpole", "fuzzy"), TRUE)

    ## coerce to list if necessary
    series <- any2list(series)

    check_consistency(series, "vltslist")

    if (!is.function(score.clus) && !(is.list(score.clus) && all(sapply(score.clus, is.function))))
        stop("Invalid evaluation function(s)")
    else if (is.list(score.clus) && !all(types %in% names(score.clus)))
        stop("The names of the 'score.clus' argument do not correspond to the provided 'types'")

    if (!is.function(pick.clus)) stop("Invalid pick function")

    ## ---------------------------------------------------------------------------------------------
    ## Some parameters
    ## ---------------------------------------------------------------------------------------------

    packages <- c("dtwclust", packages)

    dots <- list(...)

    configs <- configs[types]

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
                preproc_fun <- get0(preproc_df$preproc[i])

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
    ## Clusterings
    ## =============================================================================================

    types_objs <- mapply(configs, names(configs), seeds, SIMPLIFY = FALSE, FUN = function(config, type, seeds) {
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

        included_centroids <- c("mean", "median", "shape", "dba", "pam", "fcm", "fcmdd")

        export <- c("dots", "trace",
                    "check_consistency", "enlist", "subset_dots", "get_config_args",
                    setdiff(unique(config$preproc), "none"),
                    setdiff(unique(config$centroid), c("default", included_centroids)))

        ## -----------------------------------------------------------------------------------------
        ## perform clusterings
        ## -----------------------------------------------------------------------------------------

        config_split <- split_parallel(config, 1L)
        seeds_split <- split_parallel(seeds)
        cfg <- NULL # CHECK complains with NOTE regarding undefined global

        objs <- foreach(cfg = config_split,
                        seed = seeds_split,
                        .combine = c,
                        .multicombine = TRUE,
                        .packages = packages,
                        .export = export) %op%
                        {
                            chunk <- lapply(seq_len(nrow(cfg)), function(i) {
                                message("-------------- Using configuration: --------------")
                                print(cfg[i, , drop = FALSE])

                                ## -----------------------------------------------------------------
                                ## obtain args from configuration
                                ## -----------------------------------------------------------------

                                ## see end of this file
                                args <- get_config_args(cfg, i)

                                ## -----------------------------------------------------------------
                                ## controls for this configuration
                                ## -----------------------------------------------------------------

                                control_fun <- match.fun(paste0(type, "_control"))

                                control_args <- subset_dots(as.list(cfg[i, , drop = FALSE]),
                                                            control_fun)

                                control_args <- lapply(control_args, unlist, recursive = FALSE)

                                control <- do.call(control_fun, control_args)

                                ## -----------------------------------------------------------------
                                ## get processed series
                                ## -----------------------------------------------------------------

                                preproc_char <- cfg$preproc[i]

                                this_preproc_config <- c(list(preproc = preproc_char),
                                                         args$preproc)

                                for (this_series in series) {
                                    this_series_config <- attr(this_series, "config")
                                    if (identical(this_series_config[names(this_preproc_config)],
                                                  this_preproc_config))
                                    {
                                        break
                                    }
                                }

                                ## -----------------------------------------------------------------
                                ## distance entry to re-register in parallel worker
                                ## -----------------------------------------------------------------

                                if (type != "tadpole") {
                                    distance <- cfg$distance[i]
                                    dist_entry <- dist_entries[[distance]]

                                    if (!check_consistency(dist_entry$names[1L], "dist"))
                                        do.call(proxy::pr_DB$set_entry, dist_entry)

                                } else
                                    distance <- "dtw_basic" ## dummy

                                ## -----------------------------------------------------------------
                                ## centroid for this configuration
                                ## -----------------------------------------------------------------

                                centroid_char <- cfg$centroid[i]

                                ## -----------------------------------------------------------------
                                ## call tsclust
                                ## -----------------------------------------------------------------

                                if (centroid_char == "default") {
                                    ## do not specify centroid
                                    tsc <- do.call(tsclust,
                                                   enlist(series = this_series,
                                                          type = type,
                                                          k = cfg$k[[i]],
                                                          distance = distance,
                                                          seed = seed[[i]],
                                                          trace = trace,
                                                          args = args,
                                                          control = control,
                                                          dots = dots),
                                                   quote = TRUE)

                                } else if (centroid_char %in% included_centroids) {
                                    ## with included centroid
                                    tsc <- do.call(tsclust,
                                                   enlist(series = this_series,
                                                          type = type,
                                                          k = cfg$k[[i]],
                                                          distance = distance,
                                                          centroid = centroid_char,
                                                          seed = seed[[i]],
                                                          trace = trace,
                                                          args = args,
                                                          control = control,
                                                          dots = dots),
                                                   quote = TRUE)

                                } else {
                                    ## with centroid function
                                    delayedAssign("centroid", get0(centroid_char))

                                    tsc <- do.call(tsclust,
                                                   enlist(series = this_series,
                                                          type = type,
                                                          k = cfg$k[[i]],
                                                          distance = distance,
                                                          centroid = centroid,
                                                          seed = seed[[i]],
                                                          trace = trace,
                                                          args = args,
                                                          control = control,
                                                          dots = dots),
                                                   quote = TRUE)
                                }

                                if (inherits(tsc, "TSClusters")) {
                                    tsc@preproc <- preproc_char

                                    if (preproc_char != "none")
                                        tsc@family@preproc <- get0(preproc_char)
                                    if (centroid_char != "default")
                                        tsc@centroid <- centroid_char

                                    tsc <- list(tsc)

                                } else {
                                    tsc <- lapply(tsc, function(tsc) {
                                        tsc@preproc <- preproc_char

                                        if (preproc_char != "none")
                                            tsc@family@preproc <- get0(preproc_char)
                                        if (centroid_char != "default")
                                            tsc@centroid <- centroid_char

                                        tsc
                                    })
                                }

                                tsc
                            })

                            ## return chunk
                            unlist(chunk, recursive = FALSE, use.names = TRUE)
                        }

        ## return TSClusters
        objs
    })

    ## =============================================================================================
    ## Evaluations
    ## =============================================================================================

    if (is.function(score.clus))
        scores <- try(lapply(types_objs, score.clus, ...), silent = TRUE)
    else
        scores <- try(mapply(types_objs, score.clus[names(types_objs)],
                             SIMPLIFY = FALSE,
                             MoreArgs = dots,
                             FUN = function(objs, score_fun, ...) { score_fun(objs, ...) }),
                      silent = TRUE)

    if (inherits(scores, "try-error")) {
        warning("The score.clus function did not execute successfully.")
        scores <- NULL
        pick <- NULL

    } else {
        pick <- try(pick.clus(scores, types_objs, ...), silent = FALSE)

        if (inherits(pick, "try-error")) {
            warning("The pick.clus function did not execute successfully.")
            pick <- NULL
        }
    }

    results <- list(scores = scores, pick = pick, proc_time = proc.time() - tic)
    if (return.objects) results <- c(results, objects = types_objs)

    ## return results
    results
}

## =================================================================================================
## compare_clusterings helper: get_config_args
## =================================================================================================

get_config_args <- function(config, i) {
    preproc_args <- grepl("_preproc$", names(config))
    dist_args <- grepl("_distance$", names(config))
    cent_args <- grepl("_centroid$", names(config))

    ## --------------------------------------------------------------
    ## preprocessing args
    ## --------------------------------------------------------------

    if (config$preproc[i] != "none" && any(preproc_args)) {
        preproc_args <- as.list(config[i, preproc_args, drop = FALSE])
        names(preproc_args) <- gsub("_preproc", "", names(preproc_args))
        preproc_args <- preproc_args[!sapply(preproc_args, is.na)]

    } else
        preproc_args <- list()

    ## --------------------------------------------------------------
    ## distance args
    ## --------------------------------------------------------------

    if (any(dist_args)) {
        dist_args <- as.list(config[i, dist_args, drop = FALSE])
        names(dist_args) <- gsub("_distance", "", names(dist_args))
        dist_args <- dist_args[!sapply(dist_args, is.na)]

    } else
        dist_args <- list()

    ## --------------------------------------------------------------
    ## centroid args
    ## --------------------------------------------------------------

    if (any(cent_args)) {
        cent_args <- as.list(config[i, cent_args, drop = FALSE])
        names(cent_args) <- gsub("_preproc", "", names(cent_args))
        cent_args <- cent_args[!sapply(cent_args, is.na)]

    } else
        cent_args <- list()

    ## --------------------------------------------------------------
    ## all args
    ## --------------------------------------------------------------

    tsclust_args(preproc = preproc_args,
                 dist = dist_args,
                 cent = cent_args)
}
