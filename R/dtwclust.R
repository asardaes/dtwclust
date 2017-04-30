#' Time series clustering
#'
#' This is the *old* main function to perform time series clustering. It supports partitional,
#' hierarchical, fuzzy, k-Shape and TADPole clustering. See [tsclust()] for the new interface.
#' Please note that possible updates will only be implemented in the new function.
#'
#' @export
#'
#' @param data A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced row-wise.
#' @param type What type of clustering method to use: `"partitional"`, `"hierarchical"`, `"tadpole"`
#'   or `"fuzzy"`.
#' @param k Number of desired clusters. It may be a numeric vector with different values.
#' @param method Character vector with one or more linkage methods to use in hierarchical procedures
#'   (see [stats::hclust()]) or a function that performs hierarchical clustering based on distance
#'   matrices (e.g. [cluster::diana()]). See Hierarchical section for more details.
#' @param distance A supported distance from [proxy::dist()] (see Distance section). Ignored for
#'   `type` = `"tadpole"`.
#' @param centroid Either a supported string or an appropriate function to calculate centroids when
#'   using partitional or prototypes for hierarchical/tadpole methods. See Centroid section.
#' @param preproc Function to preprocess data. Defaults to [zscore()] *only* if `centroid` `=`
#'   `"shape"`, but will be replaced by a custom function if provided. See Preprocessing section.
#' @param dc Cutoff distance for the [TADPole()] algorithm.
#' @param control Named list of parameters or `dtwclustControl` object for clustering algorithms.
#'   See [dtwclustControl]. `NULL` means defaults.
#' @param seed Random seed for reproducibility.
#' @param distmat If a cross-distance matrix is already available, it can be provided here so it's
#'   re-used. Only relevant if `centroid` = "pam" or `type` = "hierarchical". See examples.
#' @param ... Additional arguments to pass to [proxy::dist()] or a custom function (preprocessing,
#'   centroid, etc.)
#'
#' @details
#'
#' Partitional and fuzzy clustering procedures use a custom implementation. Hierarchical clustering
#' is done with [stats::hclust()]. TADPole clustering uses the [TADPole()] function. Specifying
#' `type` = `"partitional"`, `distance` = `"sbd"` and `centroid` = `"shape"` is equivalent to the
#' k-Shape algorithm (Paparrizos and Gravano 2015).
#'
#' The `data` may be a matrix, a data frame or a list. Matrices and data frames are coerced to a
#' list, both row-wise. Only lists can have series with different lengths or multiple dimensions.
#' Most of the optimizations require series to have the same length, so consider reinterpolating
#' them to save some time (see Ratanamahatana and Keogh 2004; [reinterpolate()]). No missing values
#' are allowed.
#'
#' In the case of multivariate time series, they should be provided as a list of matrices, where
#' time spans the rows of each matrix and the variables span the columns. At the moment, only `DTW`,
#' `DTW2` and `GAK` suppport such series, which means only partitional and hierarchical procedures
#' using those distances will work. You can of course create your own custom distances. All included
#' centroid functions should work with the aforementioned format, although `shape` is **not**
#' recommended. Note that the `plot` method will simply append all dimensions (columns) one after
#' the other.
#'
#' Several parameters can be adjusted with the `control` argument. See [dtwclustControl]. In the
#' following sections, elements marked with an asterisk (*) are those that can be adjutsed with this
#' argument.
#'
#' @return
#'
#' An object with formal class [dtwclust-class].
#'
#' If `control@nrep > 1` and a partitional procedure is used, `length(method)` `> 1` and
#' hierarchical procedures are used, or `length(k)` `>` `1`, a list of objects is returned.
#'
#' @section Partitional Clustering:
#'
#'   Stochastic algorithm that creates a hard partition of the data into `k` clusters, where each
#'   cluster has a centroid. In case of time series clustering, the centroids are also time series.
#'
#'   The cluster centroids are first randomly initialized by selecting some of the series in the
#'   data. The distance between each series and each centroid is calculated, and the series are
#'   assigned to the cluster whose centroid is closest. The centroids are then updated by using a
#'   given rule, and the procedure is repeated until no series changes from one cluster to another,
#'   or until the maximum number of iterations* has been reached. The distance and centroid
#'   definitions can be specified through the corresponding parameters of this function. See their
#'   respective sections below.
#'
#'   Note that it is possible for a cluster to become empty, in which case a new cluster is
#'   reinitialized randomly. However, if you see that the algorithm doesn't converge or the overall
#'   distance sum increases, it could mean that the chosen value of `k` is too large, or the chosen
#'   distance measure is not able to assess similarity effectively. The random reinitialization
#'   attempts to enforce a certain number of clusters, but can result in instability in the
#'   aforementioned cases.
#'
#' @section Fuzzy Clustering:
#'
#'   This procedure is very similar to partitional clustering, except that each series no longer
#'   belongs exclusively to one cluster, but belongs to each cluster to a certain degree. For each
#'   series, the total degree of membership across clusters must sum to 1.
#'
#'   The default implementation uses the fuzzy c-means algorithm. In its definition, an objective
#'   function is to be minimized. The objective is defined in terms of a squared distance, which is
#'   usually the Euclidean (L2) distance, although the definition could be modified. The `distance`
#'   parameter of this function controls the one that is utilized. The fuzziness of the clustering
#'   can be controlled by means of the fuzziness exponent*. Bear in mind that the centroid
#'   definition of fuzzy c-means requires equal dimensions, which means that all series must have
#'   the same length. This problem can be circumvented by applying transformations to the series
#'   (see for example D'Urso and Maharaj (2009)).
#'
#'   Note that the fuzzy clustering could be transformed to a crisp one by finding the highest
#'   membership coefficient. Some of the slots of the object returned by this function assume this,
#'   so be careful with interpretation (see [dtwclust-class]).
#'
#' @section Hierarchical Clustering:
#'
#'   This is (by default) a deterministic algorithm that creates a hierarchy of groups by using
#'   different linkage methods (see [stats::hclust()]). The linkage method is controlled through the
#'   `method` parameter of this function, which can be a character vector with several methods, with
#'   the additional option "all" that uses all of the available methods in [stats::hclust()]. The
#'   distance to be used can be controlled with the `distance` parameter.
#'
#'   Optionally, `method` may be a **function** that performs the hierarchical clustering based on a
#'   distance matrix, such as the functions included in package \pkg{cluster}. The function will
#'   receive the `dist` object as first argument (see [stats::as.dist()]), followed by the elements
#'   in `...` that match the its formal arguments. The object it returns must support the
#'   [stats::as.hclust()] generic so that [stats::cutree()] can be used. See the examples.
#'
#'   The hierarchy does not imply a specific number of clusters, but one can be induced by cutting
#'   the resulting dendrogram (see [stats::cutree()]). This results in a crisp partition, and some
#'   of the slots of the returned object are calculated by cutting the dendrogram so that `k`
#'   clusters are created.
#'
#' @section TADPole Clustering:
#'
#'   TADPole clustering adopts a relatively new clustering framework and adapts it to time series
#'   clustering with DTW. Because of the way it works, it can be considered a kind of Partitioning
#'   Around Medoids (PAM). This means that the cluster centroids are always elements of the data.
#'   However, this algorithm is deterministic, depending on the value of the cutoff distance `dc`,
#'   which can be controlled with the corresponding parameter of this function.
#'
#'   The algorithm relies on the DTW lower bounds, which are only defined for time series of equal
#'   length. Additionally, it requires a window constraint* for DTW. See the Sakoe-Chiba constraint
#'   section below. Unlike the other algorithms, TADPole always uses DTW2 as distance (with a
#'   symmetric1 step pattern).
#'
#' @section Centroid Calculation:
#'
#'   In the case of partitional/fuzzy algorithms, a suitable function should calculate the cluster
#'   centroids at every iteration. In this case, the centroids are themselves time series. Fuzzy
#'   clustering uses the standard fuzzy c-means centroid by default.
#'
#'   In either case, a custom function can be provided. If one is provided, it will receive the
#'   following parameters with the shown names (examples for partitional clustering are shown in
#'   parenthesis):
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
#'   `shape`, `dba` and `pam` support series of different length. Also note that, for `shape` and
#'   `dba`, this support has a caveat: the final centroids' length will depend on the length of
#'   those series that were randomly chosen at the beginning of the clustering algorithm. For
#'   example, if the series in the dataset have a length of either 10 or 15, 2 clusters are desired,
#'   and the initial choice selects two series with length of 10, the final centroids will have this
#'   same length.
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
#'   via [proxy::pr_DB()], and extra parameters can be provided in `...`.
#'
#'   Note that you are free to create your own distance functions and register them. Optionally, you
#'   can use one of the following custom implementations (all registered with `proxy`):
#'
#'   - `"dtw"`: DTW, optionally with a Sakoe-Chiba/Slanted-band constraint*.
#'   - `"dtw2"`: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint*. Read
#'     details below.
#'   - `"dtw_basic"`: A custom version of DTW with less functionality, but slightly faster. See
#'     [dtw_basic()].
#'   - `"dtw_lb"`: DTW with L1 or L2 norm* and optionally a Sakoe-Chiba constraint*. Some
#'     computations are avoided by first estimating the distance matrix with Lemire's lower bound
#'     and then iteratively refining with DTW. See [dtw_lb()]. Not suitable for `pam.precompute`* =
#'     `TRUE`.
#'   - `"lbk"`: Keogh's lower bound with either L1 or L2 norm* for the Sakoe-Chiba constraint*.
#'   - `"lbi"`: Lemire's lower bound with either L1 or L2 norm* for the Sakoe-Chiba constraint*.
#'   - `"sbd"`: Shape-based distance. See [SBD()] for more details.
#'   - `"gak"`: Global alignment kernels. See [GAK()] for more details.
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
#'   If you create your own distance, register it with `proxy`, and it includes the ellipsis (`...`)
#'   in its definition, it will receive the following parameters*:
#'
#'   - `window.type`: Either `"none"` for a `NULL` `window.size`, or `"slantedband"` otherwise
#'   - `window.size`: The provided window size
#'   - `norm`: The provided desired norm
#'   - `...`: Any additional parameters provided in the original call's ellipsis
#'
#'   Whether your function makes use of them or not, is up to you.
#'
#'   If you know that the distance function is symmetric, and you use a hierarchical algorithm, or a
#'   partitional algorithm with PAM centroids and `pam.precompute`* = `TRUE`, some time can be saved
#'   by calculating only half the distance matrix. Therefore, consider setting the symmetric*
#'   control parameter to `TRUE` if this is the case.
#'
#' @section Sakoe-Chiba Constraint:
#'
#'   A global constraint to speed up the DTW calculations is the Sakoe-Chiba band (Sakoe and Chiba,
#'   1978). To use it, a window size* must be defined.
#'
#'   The windowing constraint uses a centered window. The function expects a value in `window.size`
#'   that represents the distance between the point considered and one of the edges of the window.
#'   Therefore, if, for example, `window.size = 10`, the warping for an observation \eqn{x_i}
#'   considers the points between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in `10(2) + 1 = 21`
#'   observations falling within the window.
#'
#'   The computations actually use a `slantedband` window, which is equivalent to the Sakoe-Chiba
#'   one if series have equal length, and stays along the diagonal of the local cost matrix if
#'   series have different length.
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
#'   will receive the data as first argument, followed by the contents of `...` that match its
#'   formal arguments.
#'
#'   It is convenient to provide this function if you're planning on using the [stats::predict()]
#'   generic.
#'
#' @section Repetitions:
#'
#'   Due to their stochastic nature, partitional clustering is usually repeated* several times with
#'   different random seeds to allow for different starting points. This function uses
#'   [rngtools::RNGseq()] to obtain different seed streams for each repetition, utilizing the `seed`
#'   parameter (if provided) to initialize it. If more than one repetition is made, the streams are
#'   returned in an attribute called `rng`.
#'
#'   Technically, you can also perform random repetitions for fuzzy clustering, although it might be
#'   difficult to evaluate the results, since they are usually evaluated relative to each other and
#'   not in an absolute way. Ideally, the groups wouldn't change too much once the algorithm
#'   converges.
#'
#'   Multiple values of `k` can also be provided to get different partitions using any `type` of
#'   clustering.
#'
#'   Repetitions are greatly optimized when PAM centroids are used and the whole distance matrix is
#'   precomputed*, since said matrix is reused for every repetition, and can be comptued in parallel
#'   (see Parallel section).
#'
#' @template parallel
#'
#' @section Parallel Computing:
#'
#'   Unless each repetition requires a few seconds, parallel computing probably isn't worth it. As
#'   such, I would only use this feature with `shape` and `DBA` centroids, or an expensive distance
#'   function like `DTW`.
#'
#'   If you register a parallel backend, the function will also try to do the calculation of the
#'   distance matrices in parallel. This should work with any function registered with
#'   [proxy::dist()] via [proxy::pr_DB()] whose `loop` flag is set to `TRUE`. If the function
#'   requires special packages to be loaded, provide their names in the `packages`* slot of
#'   `control`. Note that "dtwclust" is always loaded in each parallel worker, so that doesn't need
#'   to be included. Alternatively, you may want to pre-load `dtwclust` in each worker with
#'   [parallel::clusterEvalQ()].
#'
#'   In case of multiple repetitions, each worker gets a repetition task. Otherwise, the tasks
#'   (which can be a distance matrix or a centroid calculation) are usually divided into chunks
#'   according to the number of workers available.
#'
#' @section Notes:
#'
#'   The lower bounds are defined only for time series of equal length. `DTW` and `DTW2` don't
#'   require this, but they are much slower to compute.
#'
#'   The lower bounds are **not** symmetric, and `DTW` is not symmetric in general.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references.
#'
#' @seealso
#'
#' [dtwclust-methods], [dtwclust-class], [dtwclustControl], [dtwclustFamily].
#'
#' @example inst/dtwclust-examples.R
#'
dtwclust <- function(data = NULL, type = "partitional", k = 2L, method = "average",
                     distance = "dtw_basic", centroid = "pam", preproc = NULL,
                     dc = NULL, control = NULL, seed = NULL, distmat = NULL,
                     ...)
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    tic <- proc.time()

    set.seed(seed)

    if (is.null(data)) stop("No data provided")

    type <- match.arg(type, c("partitional", "hierarchical", "tadpole", "fuzzy"))

    ## coerce to list if necessary
    data <- any2list(data)

    if (any(k < 2L)) stop("At least two clusters must be defined")
    if (any(k > length(data))) stop("Cannot have more clusters than series in the data")

    MYCALL <- match.call(expand.dots = TRUE)

    ## ---------------------------------------------------------------------------------------------
    ## Control parameters
    ## ---------------------------------------------------------------------------------------------

    if (is.null(control) || is.list(control)) control <- as(control, "dtwclustControl")
    else if (class(control) != "dtwclustControl") stop("Invalid control argument")
    else methods::validObject(control)

    if (control@trace) message("Consider using the 'tsclust' function instead.\n")

    dots <- list(...)

    ## ---------------------------------------------------------------------------------------------
    ## Preprocess
    ## ---------------------------------------------------------------------------------------------

    if (!is.null(preproc) && is.function(preproc)) {
        if (has_dots(preproc)) {
            data <- preproc(data, ...)

        } else {
            data <- do.call(preproc,
                            enlist(data,
                                   dots = subset_dots(dots, preproc)))
        }

        preproc_char <- as.character(substitute(preproc))[1L]

    } else if (type == "partitional" && is.character(centroid) && centroid == "shape") {
        preproc <- zscore
        preproc_char <- "zscore"
        data <- zscore(data, ...)

    } else if (is.null(preproc)) {
        preproc <- function(x, ...) x
        preproc_char <- "none"

    } else stop("Invalid preprocessing")

    check_consistency(data, "vltslist")

    ## ---------------------------------------------------------------------------------------------
    ## Further options
    ## ---------------------------------------------------------------------------------------------

    if (control@save.data)
        datalist <- data
    else
        datalist <- list()

    diff_lengths <- different_lengths(data)

    check_consistency(distance, "dist", trace = control@trace, Lengths = diff_lengths, silent = FALSE)

    ## symmetric versions of dtw that I know of
    ## unconstrained and with symmetric1/symmetric2 is always symmetric, regardless of lengths
    ## constrained and same lengths with symmetric1/symmetric2 is also symmetric
    symmetric_pattern <- is.null(dots$step.pattern) ||
        identical(dots$step.pattern, symmetric1) ||
        identical(dots$step.pattern, symmetric2)

    if (tolower(distance) %in% c("dtw", "dtw2", "dtw_basic")) {
        control@symmetric <- symmetric_pattern && (is.null(control@window.size) || !diff_lengths)

        if (tolower(distance) == "dtw2") control@norm <- "L2"

    } else if (tolower(distance) %in% c("lbk", "lbi")) {
        control@symmetric <- FALSE

    } else if (tolower(distance) %in% c("sbd", "gak")) {
        control@symmetric <- TRUE
    }

    ## For parallel computation
    control@packages <- c("dtwclust", control@packages)

    if (type %in% c("partitional", "fuzzy")) {

        ## =========================================================================================
        ## Partitional or fuzzy
        ## =========================================================================================

        ## fuzzy default
        if (type == "fuzzy" && missing(centroid)) centroid <- "fcm"

        if (is.character(centroid)) {
            if (type == "fuzzy")
                centroid <- match.arg(centroid, c("fcm", "fcmdd"))
            else
                centroid <- match.arg(centroid, c("mean", "median", "shape", "dba", "pam"))

            ## replace any given distmat if centroid not "pam" or "fcmdd"
            if (!(centroid %in% c("pam", "fcmdd"))) distmat <- NULL

        } else {
            distmat <- NULL
        }

        if (diff_lengths) {
            if (type == "fuzzy" && is.character(centroid) && centroid == "fcm")
                stop("Fuzzy c-means does not support series with different length.")

            check_consistency(centroid, "cent", trace = control@trace)
        }

        ## -----------------------------------------------------------------------------------------
        ## Family creation, see initialization in dtwclust-methods.R
        ## -----------------------------------------------------------------------------------------

        family <- new("dtwclustFamily",
                      dist = distance,
                      allcent = centroid,
                      preproc = preproc,
                      distmat = distmat,
                      control = control,
                      fuzzy = isTRUE(type == "fuzzy"))

        if (!all(c("x", "cl_id", "k", "cent", "cl_old") %in% names(formals(family@allcent))))
            stop("The provided centroid function must have at least the following arguments with ",
                 "the shown names:\n\t",
                 paste(c("x", "cl_id", "k", "cent", "cl_old"), collapse = ", "))

        cent_char <- as.character(substitute(centroid))[1L]

        ## -----------------------------------------------------------------------------------------
        ## PAM precompute?
        ## -----------------------------------------------------------------------------------------

        ## for a check near the end, changed if appropriate
        distmat_provided <- FALSE

        ## precompute distance matrix?
        if (cent_char %in% c("pam", "fcmdd")) {
            ## check if distmat was not provided and should be precomputed
            if (!is.null(distmat)) {
                if ( nrow(distmat) != length(data) || ncol(distmat) != length(data) )
                    stop("Dimensions of provided cross-distance matrix don't correspond to length of provided data")

                ## distmat was provided in call
                distmat_provided <- TRUE

                if (control@trace) cat("\n\tDistance matrix provided...\n\n")

            } else if (control@pam.precompute || cent_char == "fcmdd") {
                if (tolower(distance) == "dtw_lb")
                    warning("Using dtw_lb with control@pam.precompute = TRUE is not advised.")

                if (control@trace) cat("\n\tPrecomputing distance matrix...\n\n")

                distmat <- do.call(family@dist, enlist(x = data, centroids = NULL, dots = dots))

                ## Redefine new distmat
                assign("distmat", distmat, environment(family@dist))
                assign("distmat", distmat, environment(family@allcent))

                gc(FALSE)
            }
        }

        ## -----------------------------------------------------------------------------------------
        ## Cluster
        ## -----------------------------------------------------------------------------------------

        if (length(k) == 1L && control@nrep == 1L) {
            rngtools::setRNG(rngtools::RNGseq(1L, seed = seed, simplify = TRUE))

            ## Just one repetition
            pc.list <- list(do.call(kcca.list,
                                    enlist(x = data,
                                           k = k,
                                           family = family,
                                           control = control,
                                           fuzzy = isTRUE(type == "fuzzy"),
                                           cent = cent_char,
                                           dots = dots)))

        } else {
            if ((foreach::getDoParName() != "doSEQ") && control@trace)
                message("Tracing of repetitions might not be available if done in parallel.\n")

            ## I need to re-register any custom distances in each parallel worker
            dist_entry <- proxy::pr_DB$get_entry(distance)

            export <- c("kcca.list", "check_consistency", "enlist")

            rng <- rngtools::RNGseq(length(k) * control@nrep, seed = seed, simplify = FALSE)
            rng0 <- lapply(parallel::splitIndices(length(rng), length(k)), function(i) rng[i])

            ## if %do% is used, the outer loop replaces value of k in this envir
            k0 <- k
            comb0 <- if (control@nrep > 1L) c else list

            i <- integer() # CHECK complains about non-initialization now

            pc.list <- foreach(k = k0, rng = rng0, .combine = comb0, .multicombine = TRUE,
                               .packages = control@packages, .export = export) %:%
                foreach(i = 1L:control@nrep, .combine = list, .multicombine = TRUE,
                        .packages = control@packages, .export = export) %op%
                        {
                            if (control@trace) message("Repetition ", i, " for k = ", k)

                            rngtools::setRNG(rng[[i]])

                            if (!check_consistency(dist_entry$names[1], "dist"))
                                do.call(proxy::pr_DB$set_entry, dist_entry)

                            pc <- do.call(kcca.list,
                                          enlist(x = data,
                                                 k = k,
                                                 family = family,
                                                 control = control,
                                                 fuzzy = isTRUE(type == "fuzzy"),
                                                 cent = cent_char,
                                                 dots = dots))

                            gc(FALSE)

                            pc
                        }
        }

        ## -----------------------------------------------------------------------------------------
        ## Prepare results
        ## -----------------------------------------------------------------------------------------

        ## Replace distmat with NULL so that, if the distance function is called again, it won't subset it
        assign("distmat", NULL, envir = environment(family@dist))

        ## If distmat was provided, let it be shown in the results
        if (distmat_provided) {
            if (is.null(attr(distmat, "method")))
                distance <- "unknown"
            else
                distance <- attr(distmat, "method")
        }

        ## Create objects
        RET <- lapply(pc.list, function(pc) {
            create_dtwclust(call = MYCALL,
                            control = control,
                            family = family,
                            distmat = distmat,

                            type = type,
                            method = "NA",
                            distance = distance,
                            centroid = cent_char,
                            preproc = preproc_char,

                            centroids = pc$centroids,
                            k = pc$k,
                            cluster = pc$cluster,
                            fcluster = pc$fcluster,
                            iter = pc$iter,
                            converged = pc$converged,
                            cldist = pc$cldist,
                            clusinfo = pc$clusinfo,

                            datalist = datalist,
                            dots = dots,
                            override.family = FALSE)
        })

        if (class(RET) != "dtwclust" && length(RET) == 1L) RET <- RET[[1L]]

    } else if (type == "hierarchical") {

        ## =========================================================================================
        ## Hierarchical
        ## =========================================================================================

        if (is.character(method)) {
            hclust_methods <- match.arg(method,
                                        c("ward.D", "ward.D2", "single", "complete",
                                          "average", "mcquitty", "median", "centroid",
                                          "all"),
                                        several.ok = TRUE)

            if (any(hclust_methods == "all"))
                hclust_methods <- c("ward.D", "ward.D2", "single", "complete",
                                    "average", "mcquitty", "median", "centroid")

        } else if (!is.function(method))
            stop("Argument 'method' must be either a supported character or a function.")

        if (tolower(distance) == "dtw_lb")
            warning("Using dtw_lb with hierarchical clustering is not advised.")

        ## -----------------------------------------------------------------------------------------
        ## Calculate distance matrix
        ## -----------------------------------------------------------------------------------------

        ## Take advantage of the function I defined for the partitional methods
        ## Which can do calculations in parallel if appropriate
        distfun <- ddist(distance = distance,
                         control = control,
                         distmat = NULL)

        if (!is.null(distmat)) {
            if ( nrow(distmat) != length(data) || ncol(distmat) != length(data) )
                stop("Dimensions of provided cross-distance matrix don't correspond to length of provided data")

            if (control@trace)
                cat("\n\tDistance matrix provided...\n")

            if (is.null(attr(distmat, "method")))
                distance <- "unknown"
            else
                distance <- attr(distmat, "method")

        } else {
            if (control@trace)
                cat("\n\tCalculating distance matrix...\n")

            ## single argument is to calculate whole distance matrix
            distmat <- do.call(distfun, enlist(x = data, centroids = NULL, dots = dots))
        }

        ## -----------------------------------------------------------------------------------------
        ## Cluster
        ## -----------------------------------------------------------------------------------------

        if (control@trace)
            cat("\n\tPerforming hierarchical clustering...\n\n")

        if (is.character(method)) {
            ## Using hclust
            hc <- lapply(hclust_methods, function(method) {
                stats::hclust(stats::as.dist(distmat), method)
            })

        } else {
            ## Using provided function
            if (has_dots(method)) {
                hc <- list(method(stats::as.dist(distmat), ...))

            } else {
                hc <- list(do.call(method,
                                   args = enlist(stats::as.dist(distmat),
                                                 dots = subset_dots(dots, method))))
            }

            method <- as.character(substitute(method))[1L]
        }

        ## Invalid centroid specifier provided?
        if (!missing(centroid) && !is.function(centroid))
            warning("The 'centroid' argument was provided but it wasn't a function, so it was ignored.")
        if (!is.function(centroid))
            centroid <- NA

        ## -----------------------------------------------------------------------------------------
        ## Prepare results
        ## -----------------------------------------------------------------------------------------

        RET <- lapply(k, function(k) {
            lapply(hc, function(hc) {
                ## cutree and corresponding centroids
                cluster <- stats::cutree(stats::as.hclust(hc), k)

                if (is.function(centroid)) {
                    allcent <- centroid
                    centroids <- lapply(1L:k, function(kcent) { centroid(data[cluster == kcent]) })
                    cent_char <- as.character(substitute(centroid))[1L]

                } else {
                    allcent <- function() {}

                    centroids <- sapply(1L:k, function(kcent) {
                        id_k <- cluster == kcent

                        d_sub <- distmat[id_k, id_k, drop = FALSE]

                        id_centroid <- which.min(apply(d_sub, 1L, sum))

                        which(id_k)[id_centroid]
                    })

                    centroids <- data[centroids]
                    cent_char <- "PAM (Hierarchical)"
                }

                create_dtwclust(stats::as.hclust(hc),
                                call = MYCALL,
                                control = control,
                                family = new("dtwclustFamily",
                                             dist = distfun,
                                             allcent = allcent,
                                             preproc = preproc),
                                distmat = distmat,

                                type = type,
                                method = if (!is.null(hc$method)) hc$method else method,
                                distance = distance,
                                centroid = cent_char,
                                preproc = preproc_char,

                                centroids = centroids,
                                k = as.integer(k),
                                cluster = cluster,

                                datalist = datalist,
                                dots = dots,
                                override.family = !is.function(centroid))
            })
        })

        RET <- unlist(RET, recursive = FALSE)

        if (class(RET) != "dtwclust" && length(RET) == 1L) RET <- RET[[1L]]

    } else if (type == "tadpole") {

        ## =========================================================================================
        ## TADPole
        ## =========================================================================================

        control@window.size <- check_consistency(control@window.size, "window")
        control@norm <- "L2"

        if (is.null(dc))
            stop("Please specify 'dc' for the TADPole algorithm")
        if (dc < 0)
            stop("The cutoff distance 'dc' must be positive")

        ## -----------------------------------------------------------------------------------------
        ## Cluster
        ## -----------------------------------------------------------------------------------------

        if (control@trace) cat("\nEntering TADPole...\n\n")

        ## mainly for predict generic
        distfun <- ddist("dtw_lb", control = control, distmat = NULL)

        ## Invalid centroid specifier provided?
        if (!missing(centroid) && !is.function(centroid))
            warning("The 'centroid' argument was provided but it wasn't a function, so it was ignored.")

        if (is.function(centroid))
            cent_char <- as.character(substitute(centroid))
        else
            cent_char <- "PAM (TADPole)"

        lb <- if (is.null(dots$lb)) "lbk" else dots$lb

        rng <- rngtools::RNGseq(length(k), seed = seed, simplify = FALSE)

        RET <- foreach(k = k, rng = rng,
                       .combine = list, .multicombine = TRUE,
                       .packages = "dtwclust", .export = "enlist") %op%
                       {
                           rngtools::setRNG(rng)

                           R <- TADPole(data, k = k, dc = dc,
                                        window.size = control@window.size, lb = lb)

                           if (control@trace) {
                               cat("TADPole completed, pruning percentage = ",
                                   formatC(100 - R$distCalcPercentage,
                                           digits = 3L,
                                           width = -1L,
                                           format = "fg"),
                                   "%\n\n",
                                   sep = "")
                           }

                           ## ----------------------------------------------------------------------
                           ## Prepare results
                           ## ----------------------------------------------------------------------

                           if (is.function(centroid)) {
                               allcent <- centroid
                               centroids <- lapply(1L:k, function(kcent) {
                                   centroid(data[R$cl == kcent])
                               })

                           } else {
                               allcent <- function(x) {}
                               centroids <- data[R$centroids]
                           }

                           obj <- create_dtwclust(call = MYCALL,
                                                  control = control,
                                                  family = new("dtwclustFamily",
                                                               dist = distfun,
                                                               allcent = allcent,
                                                               preproc = preproc),
                                                  distmat = NULL,

                                                  type = type,
                                                  distance = "dtw_lb",
                                                  centroid = cent_char,
                                                  preproc = preproc_char,

                                                  centroids = centroids,
                                                  k = as.integer(k),
                                                  cluster = as.integer(R$cl),

                                                  datalist = datalist,
                                                  dots = dots,
                                                  override.family = !is.function(centroid))

                           obj@distance <- "LB_Keogh+DTW2"
                           obj
                       }
    }

    ## =============================================================================================
    ## Finish
    ## =============================================================================================

    toc <- proc.time() - tic

    if (class(RET) == "dtwclust") {
        RET@proctime <- toc

    } else {
        RET <- lapply(RET, function(ret) {
            ret@proctime <- toc

            ret
        })
    }

    if (type %in% c("partitional", "fuzzy") && (control@nrep > 1L || length(k) > 1L))
        attr(RET, "rng") <- unlist(rng0, recursive = FALSE, use.names = FALSE)

    if (control@trace)
        cat("\tElapsed time is", toc["elapsed"], "seconds.\n\n")

    RET
}
