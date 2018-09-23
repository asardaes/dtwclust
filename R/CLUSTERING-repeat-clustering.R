#' Repeat a clustering configuration
#'
#' Repeat a clustering made with [compare_clusterings()] in order to obtain the [TSClusters-class]
#' object.
#'
#' @export
#'
#' @param series The same time series that were given to [compare_clusterings()].
#' @param clusterings The list returned by [compare_clusterings()].
#' @param config_id The character indicating which configuration should be re-computed. Obtained
#'   from the `clusterings`' `results`' data frames.
#' @param ... More arguments for [tsclust()] (e.g. `trace`).
#'
#' @details
#'
#' Since the purpose of [compare_clusterings()] is to test many configurations, it is desirable to
#' set its `return.objects` parameter to `FALSE` in order to save RAM. This function can then be
#' used to compute the clustering object for a specific `config_id`.
#'
#' @return A [TSClusters-class] object.
#'
#' @section Limitations:
#'
#'   If the preprocessing function is subject to randomness, the clustering will not be correctly
#'   re-created by this function, since [compare_clusterings()] applies all preprocessing before
#'   calling [tsclust()].
#'
#'   If any parameters were given to [compare_clusterings()] through its ellipsis, they should
#'   probably be given to this function too.
#'
repeat_clustering <- function(series, clusterings, config_id, ...) {
    if (is.null(clusterings$scores))
        stop("No scores found, are you sure you need this function?") # nocov

    results <- clusterings$results

    # get used seed
    id <- lapply(results, function(res) { which(res$config_id == config_id) })
    if (sum(lengths(id)) != 1L) stop("Configuration id not found in the clusterings' results.")
    clus_type <- names(id[which(lengths(id) == 1L)])
    top_id <- strsplit(config_id, "_")[[1L]][1L]
    seed <- clusterings$seeds[[clus_type]][[top_id]]

    # get control_args and remove them from args
    results <- results[[clus_type]]
    scores <- clusterings$scores[[clus_type]]
    score_names <- if (is.null(dim(scores))) "score" else colnames(scores)
    id <- unlist(id)
    args <- results[id, , drop = FALSE]
    keep_args <- setdiff(names(args), c("config_id", "rep", score_names))
    args <- as.list(args[keep_args])
    control_args <- names(formals(paste0(clus_type, "_control")))
    control_args <- args[names(args) %in% control_args]
    args <- args[!(names(args) %in% names(control_args))]

    # get sub_seed for non-hierarchical cases
    if (clus_type != "hierarchical") {
        handle_rngkind() # UTILS-rng.R
        matching_configs <- which(grepl(paste0(top_id, "_"), results$config_id))
        if (length(matching_configs) == 0L) matching_configs <- id
        sub_id <- which(matching_configs == id)
        .rng_ <- rng_seq(length(matching_configs), seed, simplify = FALSE)[sub_id] # UTILS-rng.R
        if (clus_type != "tadpole") .rng_ <- .rng_[[1L]]
        args$.rng_ <- .rng_
    }

    # get preproc args and remove them from args
    preproc_args <- list()
    which_preproc <- grepl("_preproc$", names(args))
    if (any(which_preproc)) {
        preproc_args <- args[which_preproc]
        names(preproc_args) <- sub("_preproc$", "", names(preproc_args))
        preproc_args <- preproc_args[!sapply(preproc_args, is.na)]
        args <- args[!which_preproc]
    }

    # get distance args and remove them from args
    distance_args <- list()
    which_distance <- grepl("_distance$", names(args))
    if (any(which_distance)) {
        distance_args <- args[which_distance]
        names(distance_args) <- sub("_distance$", "", names(distance_args))
        distance_args <- distance_args[!sapply(distance_args, is.na)]
        args <- args[!which_distance]
    }

    # get centroid args and remove them from args
    centroid_args <- list()
    which_centroid <- grepl("_centroid$", names(args))
    if (any(which_centroid)) {
        centroid_args <- args[which_centroid]
        names(centroid_args) <- sub("_centroid$", "", names(centroid_args))
        centroid_args <- centroid_args[!sapply(centroid_args, is.na)]
        args <- args[!which_centroid]
    }

    # set control args
    args$control <- do.call(paste0(clus_type, "_control"), control_args, TRUE)

    # set remaining tsclust args
    centroid_char <- args$centroid
    if (centroid_char != "default") {
        if (clus_type %in% c("hierarchical", "tadpole") || !(centroid_char %in% centroids_included))
            args$centroid <- match.fun(args$centroid)
    }
    else {
        args$centroid <- NULL
    }
    preproc_char <- if (is.null(args$preproc)) "none" else args$preproc
    args$preproc <- if (preproc_char == "none") NULL else match.fun(args$preproc)
    args$series <- series
    args$type <- clus_type
    args$seed <- seed
    args$args <- do.call(tsclust_args, quote = TRUE, args = list(
        preproc = preproc_args,
        dist = distance_args,
        cent = centroid_args)
    )

    # create TSClusters
    ret <- do.call(tsclust, c(args, list(...)), TRUE)
    ret@args <- lapply(ret@args, function(arg) { arg$.rng_ <- NULL; arg })
    ret@dots$.rng_ <- NULL
    ret@preproc <- preproc_char
    if (centroid_char != "default") ret@centroid <- centroid_char
    ret
}
