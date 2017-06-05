#' Time series clustering
#'
#' This is the new main function to perform time series clustering. It should provide the same
#' functionality as [dtwclust()], but it is hopefully more coherent in general. See the details and
#' the examples for more information, as well as the included package vignette (which can be loaded
#' by typing `vignette("dtwclust")`). A convenience wrapper is available in [compare_clusterings()].
#'
#' @export
#'
#' @param series A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced to a list row-wise.
#' @param type What type of clustering method to use: `"partitional"`, `"hierarchical"`, `"tadpole"`
#'   or `"fuzzy"`.
#' @param k Number of desired clusters. It may be a numeric vector with different values.
#' @param ... Arguments to pass to preprocessing, centroid **and** distance functions (added to
#'   `args`). Also passed to `method` from [hierarchical_control()] if it happens to be a function,
#'   and to [stats::hclust()] if it contains the `members` parameter.
#' @param preproc Function to preprocess data. Defaults to [zscore()] *only* if `centroid` `=`
#'   `"shape"`, but will be replaced by a custom function if provided.
#' @param distance A supported distance from [proxy::dist()]. Ignored for `type` `=` `"tadpole"`.
#' @param centroid Either a supported string or an appropriate function to calculate centroids when
#'   using partitional or prototypes for hierarchical/tadpole methods.
#' @param control An appropriate list of controls. See [tsclust-controls]
#' @param args An appropriate list of arguments for preprocessing, distance and centroid functions.
#'   See [tsclust_args()]
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
#' The `series` may be porovided as a matrix, a data frame or a list. Matrices and data frames are
#' coerced to a list, both row-wise. Only lists can have series with different lengths or multiple
#' dimensions. Most of the optimizations require series to have the same length, so consider
#' reinterpolating them to save some time (see Ratanamahatana and Keogh 2004; [reinterpolate()]). No
#' missing values are allowed.
#'
#' In the case of multivariate time series, they should be provided as a list of matrices, where
#' time spans the rows of each matrix and the variables span the columns (see [CharTrajMV] for an
#' example). At the moment, only `DTW`, `DTW2` and `GAK` suppport such series, which means only
#' partitional and hierarchical procedures using those distances will work. You can of course create
#' your own custom distances. All included centroid functions should work with the aforementioned
#' format, although `shape` is **not** recommended. Note that the `plot` method will simply append
#' all dimensions (columns) one after the other.
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
#'   centroids at every iteration. In this case, the centroids are themselves time series. Fuzzy
#'   clustering uses the standard fuzzy c-means centroid by default.
#'
#'   In either case, a custom function can be provided. If one is provided, it will receive the
#'   following parameters with the shown names (examples for partitional clustering are shown in
#'   parentheses):
#'
#'   - `x`: The *whole* data list (`list(ts1, ts2, ts3)`)
#'   - `cl_id`: A numeric vector with length equal to the number of series in `data`, indicating
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
#'   - "pam": Partition around medoids (PAM). This basically means that the cluster centroids are
#'     always one of the time series in the data. In this case, the distance matrix can be
#'     pre-computed once using all time series in the data and then re-used at each iteration. It
#'     usually saves overhead overall for small datasets.
#'   - "fcm": Fuzzy c-means. Only supported for fuzzy clustering and used by default in that case.
#'   - "fcmdd": Fuzzy c-medoids. Only supported for fuzzy clustering. It **always** precomputes the
#'     whole cross-distance matrix.
#'
#'   These check for the special cases where parallelization might be desired. Note that only
#'   `shape`, `dba`, `pam` and `fcmdd` support series of different length. Also note that, for
#'   `shape` and `dba`, this support has a caveat: the final centroids' length will depend on the
#'   length of those series that were randomly chosen at the beginning of the clustering algorithm.
#'   For example, if the series in the dataset have a length of either 10 or 15, 2 clusters are
#'   desired, and the initial choice selects two series with length of 10, the final centroids will
#'   have this same length.
#'
#'   As special cases, if hierarchical or tadpole clustering is used, you can provide a centroid
#'   function that takes a list of series as only input and returns a single centroid series. These
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
#'   - `"dtw"`: DTW, optionally with a Sakoe-Chiba/Slanted-band constraint.
#'   - `"dtw2"`: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint. Read
#'     details below.
#'   - `"dtw_basic"`: A custom version of DTW with less functionality, but slightly faster. See
#'     [dtw_basic()].
#'   - `"dtw_lb"`: DTW with L1 or L2 norm and optionally a Sakoe-Chiba constraint. Some
#'     computations are avoided by first estimating the distance matrix with Lemire's lower bound
#'     and then iteratively refining with DTW. See [dtw_lb()]. Not suitable for `pam.precompute` =
#'     `TRUE`.
#'   - `"lbk"`: Keogh's lower bound with either L1 or L2 norm for the Sakoe-Chiba constraint.
#'   - `"lbi"`: Lemire's lower bound with either L1 or L2 norm for the Sakoe-Chiba constraint.
#'   - `"sbd"`: Shape-based distance. See [SBD()] for more details.
#'   - `"gak"`: Global alignment kernels. See [GAK()] for more details..
#'
#'   DTW2 is done with [dtw::dtw()], but it differs from the result you would obtain if you specify
#'   `L2` as `dist.method`: with `DTW2`, pointwise distances (the local cost matrix) are calculated
#'   with `L1` norm, *each* element of the matrix is squared and the result is fed into
#'   [dtw::dtw()], which finds the optimum warping path. The square root of the resulting distance
#'   is *then* computed. See [dtw2()].
#'
#'   Only `dtw`, `dtw2`, `sbd` and `gak` support series of different length. The lower bounds are
#'   probably unsuitable for direct clustering unless series are very easily distinguishable.
#'
#'   If you know that the distance function is symmetric, and you use a hierarchical algorithm, or a
#'   partitional algorithm with PAM centroids, some time can be saved by calculating only half the
#'   distance matrix. Therefore, consider setting the symmetric control parameter to `TRUE` if this
#'   is the case.
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
#'   precomputed, since said matrix is reused for every repetition, and can be comptued in parallel
#'   (see Parallel section).
#'
#' @template parallel
#'
#' @section Parallel Computing:
#'
#'   Unless each repetition requires a few seconds, parallel computing probably isn't worth it. As
#'   such, I would only use this feature with `shape` and `DBA` centroids, or an expensive distance
#'   function like `DTW` or `GAK`.
#'
#'   If you register a parallel backend, the function will also try to do the calculation of the
#'   distance matrices in parallel. This should work with any function registered with
#'   [proxy::dist()] via [proxy::pr_DB()] whose `loop` flag is set to `TRUE`. If the function
#'   requires special packages to be loaded, provide their names in the `packages` slot of
#'   `control`. Note that "dtwclust" is always loaded in each parallel worker, so that doesn't need
#'   to be included. Alternatively, you may want to pre-load \pkg{dtwclust} in each worker with
#'   [parallel::clusterEvalQ()].
#'
#'   In case of multiple repetitions, each worker gets a repetition task. Otherwise, the tasks
#'   (which can be a distance matrix or a centroid calculation) are usually divided into chunks
#'   according to the number of workers available.
#'
#' @note
#'
#' The lower bounds are defined only for time series of equal length.
#'
#' The lower bounds are **not** symmetric, and `DTW` is not symmetric in general.
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
#' @example inst/tsclust-examples.R
#'
tsclust <- function(series = NULL, type = "partitional", k = 2L, ...,
                    preproc = NULL, distance = "dtw_basic",
                    centroid = ifelse(type == "fuzzy", "fcm", "pam"),
                    control = do.call(paste0(type, "_control"), list()),
                    args = tsclust_args(),
                    seed = NULL, trace = FALSE, error.check = TRUE)
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    tic <- proc.time()
    set.seed(seed)
    type <- match.arg(type, c("partitional", "hierarchical", "tadpole", "fuzzy"))
    series <- any2list(series) ## coerce to list if necessary

    if (any(k < 2L)) stop("At least two clusters must be defined")
    if (any(k > length(series))) stop("Cannot have more clusters than series in the dataset")
    if (!is.list(control)) stop("Invalid control argument")

    MYCALL <- match.call(expand.dots = TRUE)
    dots <- list(...)
    args <- adjust_args(args, dots) ## utils.R

    ## ---------------------------------------------------------------------------------------------
    ## Preprocess
    ## ---------------------------------------------------------------------------------------------

    if (!is.null(preproc) && is.function(preproc)) {
        series <- do.call(preproc, enlist(series, dots = subset_dots(args$preproc, preproc)))
        preproc_char <- as.character(substitute(preproc))[1L]

    } else if (type == "partitional" && is.character(centroid) && centroid == "shape") {
        preproc <- zscore
        preproc_char <- "zscore"
        series <- do.call(zscore, enlist(series, dots = args$preproc))

    } else if (is.null(preproc)) {
        preproc <- function(x, ...) { x }
        environment(preproc) <- .GlobalEnv
        preproc_char <- "none"

    } else stop("Invalid preprocessing")

    if (error.check) check_consistency(series, "vltslist")

    ## ---------------------------------------------------------------------------------------------
    ## Further options
    ## ---------------------------------------------------------------------------------------------

    ## after preprocessing!
    diff_lengths <- different_lengths(series)
    check_consistency(distance, "dist", trace = trace, diff_lengths = diff_lengths, silent = FALSE)
    cent_missing <- missing(centroid)
    cent_char <- check_consistency(centroid, "cent", clus_type = type,
                                   diff_lengths = diff_lengths, cent_missing = cent_missing)

    if (type != "tadpole") {
        ## symmetric versions of dtw that I know of
        ## unconstrained and with symmetric1/symmetric2 is always symmetric, regardless of lengths
        ## constrained and same lengths with symmetric1/symmetric2 is also symmetric
        symmetric_pattern <- is.null(args$dist$step.pattern) ||
            identical(args$dist$step.pattern, symmetric1) ||
            identical(args$dist$step.pattern, symmetric2)

        if (tolower(distance) %in% c("dtw", "dtw2", "dtw_basic"))
            control$symmetric <- symmetric_pattern && (is.null(args$dist$window.size) || !diff_lengths)
        else if (tolower(distance) %in% c("lbk", "lbi"))
            control$symmetric <- FALSE
        else if (tolower(distance) %in% c("sbd", "gak"))
            control$symmetric <- TRUE
    }

    ## pre-allocate matrices for known distances
    matrices_allocated <- FALSE
    if (type != "tadpole") {
        if (tolower(distance) %in% c("dtw_basic", "dtw_lb") && is.null(args$dist$gcm)) {
            N <- max(sapply(series, NROW))
            args$dist$gcm <- matrix(0, 2L, N + 1L)
            matrices_allocated <- TRUE

        } else if (tolower(distance) == "gak" && is.null(args$dist$logs)) {
            N <- max(sapply(series, NROW))
            args$dist$logs <- matrix(0, N + 1L, 3L)
            matrices_allocated <- TRUE
        }
    }

    ## pre-allocate matrix for DBA
    if (grepl("^dba$", cent_char, ignore.case = TRUE) && is.null(args$cent$gcm)) {
        dba_allocated <- TRUE
        if (!exists("N", mode = "integer", inherits = FALSE)) N <- max(sapply(series, NROW))
        args$cent$gcm <- matrix(0, N + 1L, N + 1L)

    } else dba_allocated <- FALSE

    RET <- switch(
        type,
        partitional =, fuzzy = {
            ## =====================================================================================
            ## Partitional or fuzzy
            ## =====================================================================================

            if (!inherits(control, "PtCtrl") && !inherits(control, "FzCtrl"))
                stop("Invalid control provided")

            nrep <- if (is.null(control$nrep)) 1L else control$nrep

            if (!is.character(centroid) || !(cent_char %in% c("pam", "fcmdd")))
                control$distmat <- NULL

            ## -------------------------------------------------------------------------------------
            ## Family creation, see initialization in tsclusters-methods.R
            ## -------------------------------------------------------------------------------------

            family <- new("tsclustFamily",
                          dist = distance,
                          allcent = centroid,
                          preproc = preproc,
                          control = control,
                          fuzzy = isTRUE(type == "fuzzy"))

            if (!all(c("x", "cl_id", "k", "cent", "cl_old") %in% names(formals(family@allcent))))
                stop("The provided centroid function must have at least the following ",
                     "arguments with the shown names:\n\t",
                     paste(c("x", "cl_id", "k", "cent", "cl_old"), collapse = ", "))

            ## -------------------------------------------------------------------------------------
            ## PAM precompute?
            ## -------------------------------------------------------------------------------------

            ## for a check near the end, changed if appropriate
            distmat_provided <- FALSE
            distmat <- control$distmat

            ## precompute distance matrix?
            if (cent_char %in% c("pam", "fcmdd")) {
                if (!is.null(distmat)) {
                    if (nrow(distmat) != length(series) || ncol(distmat) != length(series))
                        stop("Dimensions of provided cross-distance matrix don't correspond ",
                             "to length of provided data")

                    ## distmat was provided in call
                    distmat_provided <- TRUE
                    distmat <- control$distmat

                    ## see Distmat.R
                    if (!inherits(distmat, "Distmat"))
                        distmat <- Distmat$new(distmat = distmat)

                    if (trace) cat("\n\tDistance matrix provided...\n\n")

                } else if (isTRUE(control$pam.precompute) || cent_char == "fcmdd") {
                    if (tolower(distance) == "dtw_lb")
                        warning("Using dtw_lb with control$pam.precompute = TRUE is not ",
                                "advised.")

                    if (trace) cat("\n\tPrecomputing distance matrix...\n\n")

                    ## see Distmat.R
                    distmat <- Distmat$new(distmat = do.call(
                        family@dist,
                        enlist(x = series,
                               centroids = NULL,
                               dots = args$dist))
                    )

                } else {
                    if (isTRUE(control$pam.sparse) && tolower(distance) != "dtw_lb") {
                        ## see SparseDistmat.R
                        distmat <- SparseDistmat$new(series = series,
                                                     distance = distance,
                                                     control = control,
                                                     dist_args = args$dist,
                                                     error.check = FALSE)

                    } else {
                        ## see Distmat.R
                        distmat <- Distmat$new(series = series,
                                               distance = distance,
                                               control = control,
                                               dist_args = args$dist,
                                               error.check = FALSE)
                    }
                }

                ## Redefine new distmat
                control$distmat <- distmat
                environment(family@dist)$control$distmat <- distmat
                environment(family@allcent)$control$distmat <- distmat
            }

            ## -------------------------------------------------------------------------------------
            ## Cluster
            ## -------------------------------------------------------------------------------------

            if (length(k) == 1L && nrep == 1L) {
                rngtools::setRNG(rngtools::RNGseq(1L, seed = seed, simplify = TRUE))

                ## Just one repetition
                pc_list <- list(do.call(pfclust,
                                        enlist(x = series,
                                               k = k,
                                               family = family,
                                               control = control,
                                               fuzzy = isTRUE(type == "fuzzy"),
                                               cent = cent_char,
                                               trace = trace,
                                               args = args)))

            } else {
                if ((foreach::getDoParName() != "doSEQ") && trace)
                    message("Tracing of repetitions might not be available if done in ",
                            "parallel.\n")

                ## I need to re-register any custom distances in each parallel worker
                dist_entry <- proxy::pr_DB$get_entry(distance)
                export <- c("pfclust", "check_consistency", "enlist")
                rng <- rngtools::RNGseq(length(k) * nrep, seed = seed, simplify = FALSE)
                ## if %do% is used, the outer loop replaces values in this envir
                rng0 <- lapply(parallel::splitIndices(length(rng), length(k)), function(i) rng[i])
                k0 <- k
                ## sequential allows the matrix to be updated iteratively
                `%this_op%` <- if(inherits(control$distmat, "SparseDistmat")) `%do%` else `%op%`
                i <- integer() # CHECK complains about non-initialization now

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
                                    do.call(proxy::pr_DB$set_entry, dist_entry)

                                ## return
                                list(
                                    do.call(pfclust,
                                            enlist(x = series,
                                                   k = k,
                                                   family = family,
                                                   control = control,
                                                   fuzzy = isTRUE(type == "fuzzy"),
                                                   cent = cent_char,
                                                   trace = trace,
                                                   args = args))
                                )
                            }
            }

            ## -------------------------------------------------------------------------------------
            ## Prepare results
            ## -------------------------------------------------------------------------------------

            ## Replace distmat with NULL so that, if the distance function is called again,
            ## it won't subset it
            environment(family@dist)$control$distmat <- NULL

            ## If distmat was provided, let it be shown in the results
            if (distmat_provided) {
                dist_method <- attr(distmat, "method")
                distance <- if (is.null(dist_method)) "unknown" else dist_method
            }

            if (inherits(distmat, "Distmat")) distmat <- distmat$distmat
            if (matrices_allocated) { args$dist$gcm <- args$dist$logs <- NULL }
            if (dba_allocated) args$cent$gcm <- NULL

            ## Create objects
            RET <- lapply(pc_list, function(pc) {
                if (type == "partitional") {
                    new("PartitionalTSClusters",
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

                } else {
                    new("FuzzyTSClusters",
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

            ## return partitional/fuzzy
            RET
        },

        hierarchical = {
            ## =====================================================================================
            ## Hierarchical
            ## =====================================================================================

            if (!inherits(control, "HcCtrl")) stop("Invalid control provided")

            method <- control$method
            distmat <- control$distmat

            if (!is.function(centroid)) centroid <- NA
            if (tolower(distance) == "dtw_lb")
                warning("Using dtw_lb with hierarchical clustering is not advised.")

            ## -------------------------------------------------------------------------------------
            ## Calculate distance matrix
            ## -------------------------------------------------------------------------------------

            ## Take advantage of the function I defined for the partitional methods
            ## Which can do calculations in parallel if appropriate
            distfun <- ddist2(distance = distance, control = control)

            if (!is.null(distmat)) {
                if (nrow(distmat) != length(series) || ncol(distmat) != length(series))
                    stop("Dimensions of provided cross-distance matrix don't correspond to ",
                         "length of provided data")

                if (trace) cat("\n\tDistance matrix provided...\n")

                if (is.null(attr(distmat, "method")))
                    distance <- "unknown"
                else
                    distance <- attr(distmat, "method")

            } else {
                if (trace) cat("\n\tCalculating distance matrix...\n")
                distmat <- do.call(distfun, enlist(x = series, centroids = NULL, dots = args$dist))
            }

            ## -------------------------------------------------------------------------------------
            ## Cluster
            ## -------------------------------------------------------------------------------------

            if (trace) cat("\n\tPerforming hierarchical clustering...\n\n")

            if (is.character(method)) {
                ## Using hclust
                hc <- lapply(method, function(method) {
                    stats::hclust(stats::as.dist(distmat), method, members = dots$members)
                })

            } else {
                ## Using provided function
                hc <- list(do.call(method,
                                   args = enlist(stats::as.dist(distmat),
                                                 dots = subset_dots(dots, method))))

                method <- attr(method, "name")
            }

            ## -------------------------------------------------------------------------------------
            ## Prepare results
            ## -------------------------------------------------------------------------------------

            if (matrices_allocated) { args$dist$gcm <- args$dist$logs <- NULL }
            if (dba_allocated) args$cent$gcm <- NULL

            RET <- lapply(k, function(k) {
                lapply(hc, function(hc) {
                    ## cutree and corresponding centroids
                    cluster <- stats::cutree(stats::as.hclust(hc), k)

                    if (is.function(centroid)) {
                        allcent <- function(...) { list(centroid(...)) }
                        environment(allcent) <- new.env(parent = .GlobalEnv)
                        assign("centroid", centroid, environment(allcent))

                        centroids <- lapply(1L:k, function(kcent) {
                            do.call(centroid,
                                    enlist(series[cluster == kcent],
                                           dots = subset_dots(args$cent, centroid)))
                        })

                    } else {
                        allcent <- function(...) {} ## dummy

                        centroids <- sapply(1L:k, function(kcent) {
                            id_k <- cluster == kcent
                            d_sub <- distmat[id_k, id_k, drop = FALSE]
                            id_centroid <- which.min(apply(d_sub, 1L, sum))
                            which(id_k)[id_centroid]
                        })

                        centroids <- series[centroids]
                    }

                    new("HierarchicalTSClusters",
                        stats::as.hclust(hc),
                        call = MYCALL,
                        family = new("tsclustFamily",
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
            ## return hierarchical
            RET
        },

        tadpole = {
            ## =====================================================================================
            ## TADPole
            ## =====================================================================================

            if (!inherits(control, "TpCtrl")) stop("Invalid control provided")
            if (!missing(distance)) warning("The distance argument is ignored for TADPole.")

            ## -------------------------------------------------------------------------------------
            ## Parameters
            ## -------------------------------------------------------------------------------------

            ## mainly for predict generic
            distfun <- ddist2("dtw_lb", control = control)
            ## for family@dist
            args$dist$window.size <- control$window.size
            args$dist$norm <- "L2"
            args$dist$window.type <- "sakoechiba"

            ## -------------------------------------------------------------------------------------
            ## Cluster
            ## -------------------------------------------------------------------------------------

            if (trace) cat("\nEntering TADPole...\n\n")

            R <- TADPole(series,
                         k = k,
                         dc = control$dc,
                         window.size = control$window.size,
                         lb = control$lb,
                         trace = trace)

            if (length(k) == 1L && length(control$dc) == 1L) R <- list(R)

            ## -------------------------------------------------------------------------------------
            ## Prepare results
            ## -------------------------------------------------------------------------------------

            ## seeds
            rng <- rngtools::RNGseq(length(k) * length(control$dc),
                                    seed = seed,
                                    simplify = FALSE)

            ## mapply/Map were causing series to be (deeply?) copied into datalist!
            RET <- lapply(1L:length(R), function(i) {
                rngtools::setRNG(rng[[i]])
                R <- R[[i]]
                k <- length(R$centroids)

                if (is.function(centroid)) {
                    allcent <- function(...) { list(centroid(...)) }
                    environment(allcent) <- new.env(parent = .GlobalEnv)
                    assign("centroid", centroid, environment(allcent))

                    centroids <- lapply(1L:k, function(kcent) {
                        centroid(series[R$cl == kcent])
                    })

                } else {
                    allcent <- function(...) {}
                    centroids <- series[R$centroids]
                }

                if (dba_allocated) args$cent$gcm <- NULL

                obj <- new("PartitionalTSClusters",
                           call = MYCALL,
                           family = new("tsclustFamily",
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

            ## return tadpole
            RET
        }
    )

    ## =============================================================================================
    ## Finish
    ## =============================================================================================

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
