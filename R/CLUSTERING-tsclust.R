# ==================================================================================================
# Helpers
# ==================================================================================================

# Get an appropriate distance matrix object for internal use with PAM/FCMdd centroids
pam_distmat <- function(series, control, distance, cent_char, family, args, trace) {
    distmat <- control$distmat
    distmat_provided <- FALSE

    if (!is.null(distmat)) {
        if (nrow(distmat) != length(series) || ncol(distmat) != length(series))
            stop("Dimensions of provided cross-distance matrix don't correspond ",
                 "to length of provided data")
        # see S4-Distmat.R
        if (!inherits(distmat, "Distmat")) distmat <- Distmat$new(distmat = distmat)
        distmat_provided <- TRUE
        if (trace) cat("\n\tDistance matrix provided...\n\n")
    }
    else if (isTRUE(control$pam.precompute) || cent_char == "fcmdd") {
        if (distance == "dtw_lb")
            warning("Using dtw_lb with control$pam.precompute = TRUE is not ",
                    "advised.")
        if (trace) cat("\n\tPrecomputing distance matrix...\n\n")
        # see S4-Distmat.R
        distmat <- Distmat$new(distmat = do.call(
            family@dist,
            enlist(x = series,
                   centroids = NULL,
                   dots = args$dist),
            TRUE
        ))
    }
    else {
        if (isTRUE(control$pam.sparse) && distance != "dtw_lb") {
            # see S4-SparseDistmat.R
            distmat <- SparseDistmat$new(series = series,
                                         distance = distance,
                                         control = control,
                                         dist_args = args$dist,
                                         error.check = FALSE)
        }
        else {
            # see S4-Distmat.R
            distmat <- Distmat$new(series = series,
                                   distance = distance,
                                   control = control,
                                   dist_args = args$dist,
                                   error.check = FALSE)
        }
    }
    # return
    list(distmat = distmat, distmat_provided = distmat_provided)
}

# ==================================================================================================
# Main function
# ==================================================================================================

#' Time series clustering
#'
#' This is the main function to perform time series clustering. See the details and the examples for
#' more information, as well as the included package vignettes (which can be found by typing
#' `browseVignettes("dtwclust")`). A convenience wrapper is available in [compare_clusterings()].
#'
#' @export
#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom parallel splitIndices
#' @importFrom proxy pr_DB
#' @importFrom rngtools RNGseq
#' @importFrom rngtools setRNG
#' @importFrom stats as.dist
#' @importFrom stats as.hclust
#' @importFrom stats cutree
#' @importFrom stats hclust
#'
#' @param series A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced to a list row-wise (see [tslist()]).
#' @param type What type of clustering method to use: `"partitional"`, `"hierarchical"`, `"tadpole"`
#'   or `"fuzzy"`.
#' @param k Number of desired clusters. It can be a numeric vector with different values.
#' @param ... Arguments to pass to preprocessing, centroid **and** distance functions (added to
#'   `args`). Also passed to `method` from [hierarchical_control()] if it happens to be a function,
#'   and to [stats::hclust()] if it contains the `members` parameter.
#' @param preproc Function to preprocess data. Defaults to [zscore()] *only* if `centroid` `=`
#'   `"shape"`, but will be replaced by a custom function if provided.
#' @param distance A registered distance from [proxy::dist()]. Ignored for `type` `=` `"tadpole"`.
#' @param centroid Either a supported string, an appropriate function to calculate centroids when
#'   using partitional/hierarchical/tadpole methods. See Centroids section.
#' @param control An appropriate list of controls. See [tsclust-controls].
#' @param args An appropriate list of arguments for preprocessing, distance and centroid functions.
#'   See [tsclust_args()] and the examples.
#' @param seed Random seed for reproducibility.
#' @param trace Logical flag. If `TRUE`, more output regarding the progress is printed to screen.
#' @template error-check
#'
#' @details
#'
#' Partitional and fuzzy clustering procedures use a custom implementation. Hierarchical clustering
#' is done with [stats::hclust()] by default. TADPole clustering uses the [TADPole()] function.
#' Specifying `type` = `"partitional"`, `preproc` = `zscore`, `distance` = `"sbd"` and `centroid` =
#' `"shape"` is equivalent to the k-Shape algorithm (Paparrizos and Gravano 2015).
#'
#' The `series` may be provided as a matrix, a data frame or a list. Matrices and data frames are
#' coerced to a list, both row-wise. Only lists can have series with different lengths or multiple
#' dimensions. Most of the optimizations require series to have the same length, so consider
#' reinterpolating them to save some time (see Ratanamahatana and Keogh 2004; [reinterpolate()]). No
#' missing values are allowed.
#'
#' In the case of multivariate time series, they should be provided as a list of matrices, where
#' time spans the rows of each matrix and the variables span the columns (see [CharTrajMV] for an
#' example). At the moment, only `dtw_basic`, `DTW`, `DTW2` and `GAK` support such series. You can
#' of course create your own custom distances. All included centroid functions should work with the
#' aforementioned format, although `shape` is *not* recommended. Note that the `plot` method will
#' simply append all dimensions (columns) one after the other.
#'
#' @return
#'
#' An object with an appropriate class from [TSClusters-class].
#'
#' If `control@nrep > 1` and a partitional procedure is used, `length(method)` `> 1` and
#' hierarchical procedures are used, or `length(k)` `>` `1`, a list of objects is returned.
#'
#' @section Centroid Calculation:
#'
#'   In the case of partitional/fuzzy algorithms, a suitable function should calculate the cluster
#'   centroids at every iteration. In this case, the centroids may also be time series. Fuzzy
#'   clustering uses the standard fuzzy c-means centroid by default.
#'
#'   In either case, a custom function can be provided. If one is provided, it will receive the
#'   following parameters with the shown names (examples for partitional clustering are shown in
#'   parentheses):
#'
#'   - `x`: The *whole* data list (`list(ts1, ts2, ts3)`)
#'   - `cl_id`: An integer vector with length equal to the number of series in `data`, indicating
#'     which cluster a series belongs to (`c(1L, 2L, 2L)`)
#'   - `k`: The desired number of total clusters (`2L`)
#'   - `cent`: The current centroids in order, in a list (`list(centroid1, centroid2)`)
#'   - `cl_old`: The membership vector of the *previous* iteration (`c(1L, 1L, 2L)`)
#'   - The elements of `...` that match its formal arguments
#'
#'   In case of fuzzy clustering, the membership vectors (2nd and 5th elements above) are matrices
#'   with number of rows equal to amount of elements in the data, and number of columns equal to the
#'   number of desired clusters. Each row must sum to 1.
#'
#'   The other option is to provide a character string for the custom implementations. The following
#'   options are available:
#'
#'   - "mean": The average along each dimension. In other words, the average of all \eqn{x^j_i}
#'   among the \eqn{j} series that belong to the same cluster for all time points \eqn{t_i}.
#'   - "median": The median along each dimension. Similar to mean.
#'   - "shape": Shape averaging. By default, all series are z-normalized in this case, since the
#'     resulting centroids will also have this normalization. See [shape_extraction()] for more
#'     details.
#'   - "dba": DTW Barycenter Averaging. See [DBA()] for more details.
#'   - "sdtw_cent": Soft-DTW centroids, See [sdtw_cent()] for more details.
#'   - "pam": Partition around medoids (PAM). This basically means that the cluster centroids are
#'     always one of the time series in the data. In this case, the distance matrix can be
#'     pre-computed once using all time series in the data and then re-used at each iteration. It
#'     usually saves overhead overall for small datasets (see [tsclust-controls]).
#'   - "fcm": Fuzzy c-means. Only supported for fuzzy clustering and used by default in that case.
#'   - "fcmdd": Fuzzy c-medoids. Only supported for fuzzy clustering. It **always** precomputes the
#'     whole cross-distance matrix.
#'
#'   The `dba`, `shape` and `sdtw_cent` implementations check for parallelization. Note that only
#'   `shape`, `dba`, `sdtw_cent`, `pam` and `fcmdd` support series of different length. Also note
#'   that for `shape`, `dba` and `sdtw_cent`, this support has a caveat: the final centroids' length
#'   will depend on the length of those series that were randomly chosen at the beginning of the
#'   clustering algorithm. For example, if the series in the dataset have a length of either 10 or
#'   15, 2 clusters are desired, and the initial choice selects two series with length of 10, the
#'   final centroids will have this same length.
#'
#'   As special cases, if hierarchical or tadpole clustering is used, you can provide a centroid
#'   function that takes a list of series as first input. It will also receive the contents of
#'   `args$cent` that match its formal arguments, and should return a single centroid series. These
#'   centroids are returned in the `centroids` slot. By default, a type of PAM centroid function is
#'   used.
#'
#' @section Distance Measures:
#'
#'   The distance measure to be used with partitional, hierarchical and fuzzy clustering can be
#'   modified with the `distance` parameter. The supported option is to provide a string, which must
#'   represent a compatible distance registered with `proxy`'s [proxy::dist()]. Registration is done
#'   via [proxy::pr_DB()], and extra parameters can be provided in `args$dist` (see the examples).
#'
#'   Note that you are free to create your own distance functions and register them. Optionally, you
#'   can use one of the following custom implementations (all registered with `proxy`):
#'
#'   - `"dtw"`: DTW, optionally with a Sakoe-Chiba/Slanted-band constraint. Done with [dtw::dtw()].
#'   - `"dtw2"`: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint. See
#'     [dtw2()].
#'   - `"dtw_basic"`: A custom version of DTW with less functionality, but faster. See
#'     [dtw_basic()].
#'   - `"dtw_lb"`: DTW with L1 or L2 norm and optionally a Sakoe-Chiba constraint. Some
#'     computations are avoided by first estimating the distance matrix with Lemire's lower bound
#'     and then iteratively refining with DTW. See [dtw_lb()]. Not suitable for `pam.precompute` =
#'     `TRUE` nor hierarchical clustering.
#'   - `"lbk"`: Keogh's lower bound for DTW with either L1 or L2 norm for the Sakoe-Chiba
#'     constraint. See [lb_keogh()].
#'   - `"lbi"`: Lemire's lower bound for DTW with either L1 or L2 norm for the Sakoe-Chiba
#'     constraint. See [lb_improved()].
#'   - `"sbd"`: Shape-based distance. See [SBD()].
#'   - `"gak"`: Global alignment kernels. See [GAK()].
#'   - `"sdtw"`: Soft-DTW. See [sdtw()].
#'
#'   Out of the aforementioned, only the distances based on DTW lower bounds *don't* support series
#'   of different length. The lower bounds are probably unsuitable for direct clustering unless
#'   series are very easily distinguishable.
#'
#'   If you know that the distance function is symmetric, and you use a hierarchical algorithm, or a
#'   partitional algorithm with PAM centroids, or fuzzy c-medoids, some time can be saved by
#'   calculating only half the distance matrix. Therefore, consider setting the symmetric control
#'   parameter to `TRUE` if this is the case (see [tsclust-controls]).
#'
#' @section Preprocessing:
#'
#'   It is strongly advised to use z-normalization in case of `centroid = "shape"`, because the
#'   resulting series have this normalization (see [shape_extraction()]). Therefore, [zscore()] is
#'   the default in this case. The user can, however, specify a custom function that performs any
#'   transformation on the data, but the user must make sure that the format stays consistent, i.e.
#'   a list of time series.
#'
#'   Setting to `NULL` means no preprocessing (except for `centroid = "shape"`). A provided function
#'   will receive the data as first argument, followed by the contents of `args$preproc` that match
#'   its formal arguments.
#'
#'   It is convenient to provide this function if you're planning on using the [stats::predict()]
#'   generic (see also [tsclusters-methods]).
#'
#' @section Repetitions:
#'
#'   Due to their stochastic nature, partitional clustering is usually repeated several times with
#'   different random seeds to allow for different starting points. This function uses
#'   [rngtools::RNGseq()] to obtain different seed streams for each repetition, utilizing the `seed`
#'   parameter (if provided) to initialize it. If more than one repetition is made, the streams are
#'   returned in an attribute called `rng`.
#'
#'   Multiple values of `k` can also be provided to get different partitions using any `type` of
#'   clustering.
#'
#'   Repetitions are greatly optimized when PAM centroids are used and the whole distance matrix is
#'   precomputed, since said matrix is reused for every repetition.
#'
#' @template parallel
#'
#' @section Parallel Computing:
#'
#'   Multi-processing is used in partitional and fuzzy clustering for multiple values of `k` and/or
#'   `nrep` (in [partitional_control()]). See [TADPole()] to know how it uses parallelization. For
#'   cross-distance matrix calculations, the parallelization strategy depends on whether the
#'   distance is included with \pkg{dtwclust} or not, see the caveats in [tsclustFamily-class].
#'
#'   If you register a parallel backend and special packages must be loaded, provide their names in
#'   the `packages` element of `control`. Note that "dtwclust" is always loaded in each parallel
#'   worker, so that doesn't need to be included. Alternatively, you may want to pre-load
#'   \pkg{dtwclust} in each worker with [parallel::clusterEvalQ()].
#'
#' @note
#'
#' The lower bounds are defined only for time series of equal length. They are **not** symmetric,
#' and `DTW` is not symmetric in general.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references (which can be loaded by typing
#' `vignette("dtwclust")`).
#'
#' @seealso
#'
#' [TSClusters-class], [tsclusters-methods], [tsclustFamily-class], [tsclust-controls],
#' [compare_clusterings()].
#'
#' @example man-examples/tsclust-examples.R
#'
tsclust <- function(series = NULL, type = "partitional", k = 2L, ...,
                    preproc = NULL, distance = "dtw_basic",
                    centroid = ifelse(type == "fuzzy", "fcm", "pam"),
                    control = do.call(paste0(type, "_control"), list()),
                    args = tsclust_args(),
                    seed = NULL, trace = FALSE, error.check = TRUE)
{
    # ==============================================================================================
    # Start
    # ==============================================================================================

    tic <- proc.time()
    set.seed(seed)
    type <- match.arg(type, c("partitional", "hierarchical", "tadpole", "fuzzy"))
    series <- tslist(series, error.check) # coerce to list if necessary
    if (any(k < 2L)) stop("At least two clusters must be defined")
    if (any(k > length(series))) stop("Cannot have more clusters than series in the dataset")
    if (!is.list(control)) stop("Invalid control argument")
    MYCALL <- match.call(expand.dots = TRUE)
    dots <- list(...)
    args <- adjust_args(args, dots) # UTILS-utils.R

    # ----------------------------------------------------------------------------------------------
    # Preprocess
    # ----------------------------------------------------------------------------------------------

    if (!is.null(preproc) && is.function(preproc)) {
        series <- do.call(preproc, enlist(series, dots = subset_dots(args$preproc, preproc)), TRUE)
        preproc_char <- as.character(substitute(preproc))[1L]
    }
    else if (type == "partitional" && is.character(centroid) && centroid == "shape") {
        preproc <- zscore
        preproc_char <- "zscore"
        series <- do.call(zscore, enlist(series, dots = args$preproc), TRUE)
    }
    else if (is.null(preproc)) {
        preproc <- function(x, ...) { x } # nocov
        environment(preproc) <- .GlobalEnv
        preproc_char <- "none"
    }
    else stop("Invalid preprocessing")

    if (error.check) check_consistency(series, "vltslist")

    # ----------------------------------------------------------------------------------------------
    # Further options
    # ----------------------------------------------------------------------------------------------

    # after preprocessing!
    distance_missing <- missing(distance)
    diff_lengths <- different_lengths(series)
    check_consistency(distance, "dist", trace = trace, diff_lengths = diff_lengths, silent = FALSE)
    distance <- tolower(distance)
    cent_missing <- missing(centroid)
    cent_char <- check_consistency(centroid, "cent", clus_type = type,
                                   diff_lengths = diff_lengths, cent_missing = cent_missing,
                                   cent_char = as.character(substitute(centroid))[1L])

    if (type != "tadpole") {
        # symmetric versions of dtw that I know of
        # unconstrained and with symmetric1/symmetric2 is always symmetric, regardless of lengths
        # constrained and same lengths with symmetric1/symmetric2 is also symmetric
        symmetric_pattern <- is.null(args$dist$step.pattern) ||
            identical(args$dist$step.pattern, dtw::symmetric1) ||
            identical(args$dist$step.pattern, dtw::symmetric2)

        if (distance %in% c("dtw", "dtw2", "dtw_basic"))
            control$symmetric <- symmetric_pattern && (is.null(args$dist$window.size) || !diff_lengths)
        else if (distance %in% c("lbk", "lbi"))
            control$symmetric <- FALSE
        else if (distance %in% c("sbd", "gak", "sdtw"))
            control$symmetric <- TRUE

        if (distance == "dtw_lb" && isTRUE(args$dist$nn.margin != 1L)) { # nocov start
            warning("Using dtw_lb in tsclust() always uses row-wise nearest neighbors.")
            args$dist$nn.margin <- 1L
        } # nocov end
    }

    RET <- switch(
        type,
        partitional =, fuzzy = {
            # ======================================================================================
            # Partitional or fuzzy
            # ======================================================================================

            if (!inherits(control, "PtCtrl") && !inherits(control, "FzCtrl"))
                stop("Invalid control provided")
            nrep <- if (is.null(control$nrep)) 1L else control$nrep
            if (!is.character(centroid) || !(cent_char %in% c("pam", "fcmdd")))
                control$distmat <- NULL

            # --------------------------------------------------------------------------------------
            # Family creation, see initialization in S4-tsclustFamily.R
            # --------------------------------------------------------------------------------------

            family <- methods::new("tsclustFamily",
                                   dist = distance,
                                   allcent = centroid,
                                   preproc = preproc,
                                   control = control,
                                   fuzzy = isTRUE(type == "fuzzy"))

            if (!all(c("x", "cl_id", "k", "cent", "cl_old") %in% names(formals(family@allcent))))
                stop("The provided centroid function must have at least the following ",
                     "arguments with the shown names:\n\t",
                     paste(c("x", "cl_id", "k", "cent", "cl_old"), collapse = ", "))

            # --------------------------------------------------------------------------------------
            # PAM precompute?
            # --------------------------------------------------------------------------------------

            # precompute distance matrix?
            if (cent_char %in% c("pam", "fcmdd")) {
                dm <- pam_distmat(series, control, distance, cent_char, family, args, trace)
                distmat <- dm$distmat
                distmat_provided <- dm$distmat_provided
                # Redefine new distmat if appropriate
                control$distmat <- distmat
                environment(family@allcent)$control$distmat <- distmat
                if (!(distance == "dtw_lb" && !isTRUE(control$pam.precompute)))
                    environment(family@dist)$control$distmat <- distmat
            }
            else {
                distmat <- NULL
                distmat_provided <- FALSE
            }

            # --------------------------------------------------------------------------------------
            # Cluster
            # --------------------------------------------------------------------------------------

            if (length(k) == 1L && nrep == 1L) {
                rngtools::setRNG(rngtools::RNGseq(1L, seed = seed, simplify = TRUE))
                # just one repetition,
                # done like this so dist/cent functions can take advantage of parallelization
                pc_list <- list(do.call(pfclust,
                                        enlist(x = series,
                                               k = k,
                                               family = family,
                                               control = control,
                                               fuzzy = isTRUE(type == "fuzzy"),
                                               cent = cent_char,
                                               trace = trace,
                                               args = args),
                                        TRUE))
            }
            else {
                # I need to re-register any custom distances in each parallel worker
                dist_entry <- proxy::pr_DB$get_entry(distance)
                export <- c("pfclust", "check_consistency", "enlist")
                rng <- rngtools::RNGseq(length(k) * nrep, seed = seed, simplify = FALSE)
                # if %do% is used, the outer loop replaces values in this envir
                rng0 <- lapply(parallel::splitIndices(length(rng), length(k)), function(i) { rng[i] })
                k0 <- k
                # sequential allows the matrix to be updated iteratively
                `%this_op%` <- if(inherits(control$distmat, "SparseDistmat")) `%do%` else `%op%`
                i <- integer() # CHECK complains about non-initialization
                # cluster
                pc_list <- foreach(k = k0, rng = rng0,
                                   .combine = c, .multicombine = TRUE,
                                   .packages = control$packages,
                                   .export = export) %:%
                    foreach(i = 1L:nrep,
                            .combine = c, .multicombine = TRUE,
                            .packages = control$packages,
                            .export = export) %this_op%
                            {
                                if (trace) message("Repetition ", i, " for k = ", k)
                                rngtools::setRNG(rng[[i]])

                                if (!check_consistency(dist_entry$names[1L], "dist"))
                                    do.call(proxy::pr_DB$set_entry, dist_entry, TRUE)

                                # return
                                list(
                                    do.call(pfclust,
                                            enlist(x = series,
                                                   k = k,
                                                   family = family,
                                                   control = control,
                                                   fuzzy = isTRUE(type == "fuzzy"),
                                                   cent = cent_char,
                                                   trace = trace,
                                                   args = args),
                                            TRUE)
                                )
                            }
            }

            # --------------------------------------------------------------------------------------
            # Prepare results
            # --------------------------------------------------------------------------------------

            # Replace distmat with NULL so that, if the distance function is called again,
            # it won't subset it
            environment(family@dist)$control$distmat <- NULL
            # If distmat was provided, let it be shown in the results
            if (distmat_provided) {
                dist_method <- attr(distmat, "method")
                distance <- if (is.null(dist_method)) "unknown" else dist_method
            }
            if (inherits(distmat, "Distmat")) distmat <- distmat$distmat
            # Create objects
            RET <- lapply(pc_list, function(pc) {
                if (type == "partitional") {
                    methods::new("PartitionalTSClusters",
                                 call = MYCALL,
                                 family = family,
                                 control = control,
                                 datalist = series,

                                 type = type,
                                 distance = distance,
                                 centroid = cent_char,
                                 preproc = preproc_char,

                                 k = pc$k,
                                 cluster = pc$cluster,
                                 centroids = pc$centroids,
                                 distmat = distmat,

                                 dots = dots,
                                 args = args,

                                 iter = pc$iter,
                                 converged = pc$converged,
                                 clusinfo = pc$clusinfo,
                                 cldist = pc$cldist,

                                 override.family = FALSE)
                }
                else {
                    methods::new("FuzzyTSClusters",
                                 call = MYCALL,
                                 family = family,
                                 control = control,
                                 datalist = series,

                                 type = type,
                                 distance = distance,
                                 centroid = cent_char,
                                 preproc = preproc_char,

                                 k = pc$k,
                                 cluster = pc$cluster,
                                 centroids = pc$centroids,
                                 distmat = distmat,

                                 dots = dots,
                                 args = args,

                                 iter = pc$iter,
                                 converged = pc$converged,
                                 fcluster = pc$fcluster,

                                 override.family = FALSE)
                }
            })
            if (any(!sapply(RET, methods::slot, name = "converged")))
                warning("At least one clustering did not converge within the allowed iterations.")
            # return partitional/fuzzy
            RET
        },
        hierarchical = {
            # ======================================================================================
            # Hierarchical
            # ======================================================================================

            if (!inherits(control, "HcCtrl")) stop("Invalid control provided")
            method <- control$method
            distmat <- control$distmat
            if (!is.function(centroid)) centroid <- NA
            if (distance == "dtw_lb")
                warning("Using dtw_lb with hierarchical clustering is not advised.")

            # --------------------------------------------------------------------------------------
            # Calculate distance matrix
            # --------------------------------------------------------------------------------------

            # Take advantage of the function I defined for the partitional methods
            # Which can do calculations in parallel if appropriate
            distfun <- ddist2(distance = distance, control = control)

            if (!is.null(distmat)) {
                if (nrow(distmat) != length(series) || ncol(distmat) != length(series))
                    stop("Dimensions of provided cross-distance matrix don't correspond to ",
                         "length of provided data")
                if (trace) cat("\nDistance matrix provided...\n")

                if (is.null(attr(distmat, "method")))
                    stop("Provided distance matrix does not include the 'method' attribute")
                else
                    distance <- attr(distmat, "method")
            }
            else {
                if (trace) cat("\nCalculating distance matrix...\n")
                distmat <- do.call(distfun,
                                   enlist(x = series, centroids = NULL, dots = args$dist),
                                   TRUE)
            }

            # --------------------------------------------------------------------------------------
            # Cluster
            # --------------------------------------------------------------------------------------

            if (trace) cat("Performing hierarchical clustering...\n")
            if (!base::isSymmetric(base::as.matrix(distmat)))
                warning("Distance matrix is not symmetric, ",
                        "and hierarchical clustering assumes it is ",
                        "(it ignores the upper triangular).")
            if (is.character(method)) {
                # Using hclust
                hc <- lapply(method, function(method) {
                    stats::hclust(stats::as.dist(distmat), method, members = dots$members)
                })
            }
            else {
                # Using provided function
                hc <- list(do.call(method,
                                   args = enlist(stats::as.dist(distmat),
                                                 dots = subset_dots(dots, method)),
                                   TRUE))
                method <- attr(method, "name")
            }

            # --------------------------------------------------------------------------------------
            # Prepare results
            # --------------------------------------------------------------------------------------

            if (trace) cat("Extracting centroids...\n\n")
            RET <- lapply(k, function(k) {
                lapply(hc, function(hc) {
                    # cutree and corresponding centroids
                    cluster <- stats::cutree(stats::as.hclust(hc), k)
                    if (is.function(centroid)) {
                        allcent <- function(...) { list(centroid(...)) }
                        environment(allcent) <- new.env(parent = .GlobalEnv)
                        assign("centroid", centroid, environment(allcent))
                        centroids <- lapply(1L:k, function(kcent) {
                            do.call(centroid,
                                    enlist(series[cluster == kcent],
                                           dots = subset_dots(args$cent, centroid)),
                                    TRUE)
                        })
                    }
                    else {
                        allcent <- function(...) {} # dummy
                        centroids <- sapply(1L:k, function(kcent) {
                            id_k <- cluster == kcent
                            d_sub <- distmat[id_k, id_k, drop = FALSE]
                            id_centroid <- which.min(apply(d_sub, 1L, sum))
                            which(id_k)[id_centroid]
                        })
                        centroids <- series[centroids]
                    }

                    methods::new("HierarchicalTSClusters",
                                 stats::as.hclust(hc),
                                 call = MYCALL,
                                 family = methods::new("tsclustFamily",
                                                       dist = distfun,
                                                       allcent = allcent,
                                                       preproc = preproc),
                                 control = control,
                                 datalist = series,

                                 type = type,
                                 distance = distance,
                                 centroid = cent_char,
                                 preproc = preproc_char,

                                 k = as.integer(k),
                                 cluster = cluster,
                                 centroids = centroids,
                                 distmat = distmat,

                                 dots = dots,
                                 args = args,

                                 method = if (!is.null(hc$method)) hc$method else method,

                                 override.family = !is.function(centroid))
                })
            })
            RET <- unlist(RET, recursive = FALSE)
            # return hierarchical
            RET
        },
        tadpole = {
            # ======================================================================================
            # TADPole
            # ======================================================================================

            if (!inherits(control, "TpCtrl")) stop("Invalid control provided")
            if (!distance_missing) warning("The distance argument is ignored for TADPole.")

            # --------------------------------------------------------------------------------------
            # Parameters
            # --------------------------------------------------------------------------------------

            # mainly for predict generic
            distfun <- ddist2("dtw_lb", control = control)
            # for family@dist
            args$dist$window.size <- control$window.size
            args$dist$norm <- "L2"
            args$dist$window.type <- "sakoechiba"

            # --------------------------------------------------------------------------------------
            # Cluster
            # --------------------------------------------------------------------------------------

            if (trace) cat("\n\tEntering TADPole...\n\n")
            R <- TADPole(series,
                         k = k,
                         dc = control$dc,
                         window.size = control$window.size,
                         lb = control$lb,
                         trace = trace)
            if (length(k) == 1L && length(control$dc) == 1L) R <- list(R)

            # --------------------------------------------------------------------------------------
            # Prepare results
            # --------------------------------------------------------------------------------------

            # seeds
            rng <- rngtools::RNGseq(length(k) * length(control$dc),
                                    seed = seed,
                                    simplify = FALSE)

            RET <- Map(R, rng, f = function(R, rng) {
                rngtools::setRNG(rng)
                k <- length(R$centroids)
                if (is.function(centroid)) {
                    allcent <- function(...) { list(centroid(...)) }
                    environment(allcent) <- new.env(parent = .GlobalEnv)
                    assign("centroid", centroid, environment(allcent))
                    centroids <- lapply(1L:k, function(kcent) {
                        do.call(centroid,
                                enlist(series[R$cl == kcent],
                                       dots = subset_dots(args$cent, centroid)),
                                TRUE)
                    })
                }
                else {
                    allcent <- function(...) {}
                    centroids <- series[R$centroids]
                }

                obj <- methods::new("PartitionalTSClusters",
                                    call = MYCALL,
                                    family = methods::new("tsclustFamily",
                                                          dist = distfun,
                                                          allcent = allcent,
                                                          preproc = preproc),
                                    control = control,
                                    datalist = series,

                                    type = type,
                                    distance = "dtw_lb",
                                    centroid = cent_char,
                                    preproc = preproc_char,

                                    k = as.integer(k),
                                    cluster = R$cl,
                                    centroids = centroids,
                                    distmat = NULL,

                                    dots = dots,
                                    args = args,

                                    override.family = !is.function(centroid))
                obj@distance <- "LB+DTW2"
                obj
            })
            # return tadpole
            RET
        }
    )

    # ==============================================================================================
    # Finish
    # ==============================================================================================

    toc <- proc.time() - tic
    RET <- lapply(RET, function(ret) {
        ret@proctime <- toc
        ret@seed <- as.integer(seed)
        ret
    })
    if (length(RET) == 1L)
        RET <- RET[[1L]]
    else if (type %in% c("partitional", "fuzzy"))
        attr(RET, "rng") <- unlist(rng0, recursive = FALSE, use.names = FALSE)
    if (trace) cat("\tElapsed time is", toc["elapsed"], "seconds.\n\n")
    RET
}
