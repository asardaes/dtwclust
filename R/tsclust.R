#' Time series clustering
#'
#' This is the new experimental main function to perform time series clustering. It should provide
#' the same functionality as \code{\link{dtwclust}}, but it is hopefully more coherent in general.
#' \strong{For now, it is subject to change}. Feedback is appreciated at
#' \url{https://github.com/asardaes/dtwclust/issues}.
#'
#' @export
#'
#' @param series A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced row-wise.
#' @param type What type of clustering method to use: \code{"partitional"}, \code{"hierarchical"},
#'   \code{"tadpole"} or \code{"fuzzy"}.
#' @param k Number of desired clusters. It may be a numeric vector with different values.
#' @param ... Arguments for other functions. Currently passed to \code{method} from
#'   \code{\link{hierarchical_control}} if it happens to be a function.
#' @param preproc Function to preprocess data. Defaults to \code{\link{zscore}} \emph{only} if
#'   \code{centroid} \code{=} \code{"shape"}, but will be replaced by a custom function if provided.
#' @param distance  A supported distance from \code{proxy}'s \code{\link[proxy]{dist}}. Ignored for
#'   \code{type} \code{=} \code{"tadpole"}.
#' @param centroid Either a supported string or an appropriate function to calculate centroids when
#'   using partitional or prototypes for hierarchical/tadpole methods.
#' @param control An appropriate list of controls. See \code{\link{tsclust-controls}}
#' @param args An appropriate list of arguments for preprocessing, distance and centroid functions.
#'   See \code{\link{tsclust_args}}
#' @param seed Random seed for reproducibility.
#' @param trace Logical flag. If \code{TRUE}, more output regarding the progress is printed to
#'   screen.
#'
#' @details
#'
#' Partitional and fuzzy clustering procedures use a custom implementation. Hierarchical clustering
#' is done with \code{\link[stats]{hclust}}. TADPole clustering uses the \code{\link{TADPole}}
#' function. Specifying \code{type} = \code{"partitional"}, \code{distance} = \code{"sbd"} and
#' \code{centroid} = \code{"shape"} is equivalent to the k-Shape algorithm (Paparrizos and Gravano
#' 2015).
#'
#' The \code{data} may be a matrix, a data frame or a list. Matrices and data frames are coerced to
#' a list, both row-wise. Only lists can have series with different lengths or multiple dimensions.
#' Most of the optimizations require series to have the same length, so consider reinterpolating
#' them to save some time (see Ratanamahatana and Keogh 2004; \code{\link{reinterpolate}}). No
#' missing values are allowed.
#'
#' In the case of multivariate time series, they should be provided as a list of matrices, where
#' time spans the rows of each matrix and the variables span the columns. At the moment, only
#' \code{DTW}, \code{DTW2} and \code{GAK} suppport such series, which means only partitional and
#' hierarchical procedures using those distances will work. You can of course create your own custom
#' distances. All included centroid functions should work with the aforementioned format, although
#' \code{shape} is \strong{not} recommended. Note that the \code{plot} method will simply append all
#' dimensions (columns) one after the other.
#'
#' @return
#'
#' An object with an appropriate class from \code{\link{TSClusters-class}}.
#'
#' If \code{control@nrep > 1} and a partitional procedure is used, \code{length(method)} \code{> 1}
#' and hierarchical procedures are used, or \code{length(k)} \code{>} \code{1}, a list of objects is
#' returned.
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
#'   \itemize{
#'     \item \code{"x"}: The \emph{whole} data list (\code{list(ts1, ts2, ts3)})
#'     \item \code{"cl_id"}: A numeric vector with length equal to the number of series in
#'       \code{data}, indicating which cluster a series belongs to (\code{c(1L, 2L, 2L)})
#'     \item \code{"k"}: The desired number of total clusters (\code{2L})
#'     \item \code{"cent"}: The current centroids in order, in a list (\code{list(centroid1,
#'       centroid2)})
#'     \item \code{"cl_old"}: The membership vector of the \emph{previous} iteration (\code{c(1L,
#'       1L, 2L)})
#'     \item The elements of \code{args$cent} that match its formal arguments
#'   }
#'
#'   In case of fuzzy clustering, the membership vectors (2nd and 5th elements above) are matrices
#'   with number of rows equal to amount of elements in the data, and number of columns equal to the
#'   number of desired clusters. Each row must sum to 1.
#'
#'   The other option is to provide a character string for the custom implementations. The following
#'   options are available:
#'
#'   \itemize{
#'     \item "mean": The average along each dimension. In other words, the average of all
#'       \eqn{x^j_i} among the \eqn{j} series that belong to the same cluster for all time points
#'       \eqn{t_i}.
#'     \item "median": The median along each dimension. Similar to mean.
#'     \item "shape": Shape averaging. By default, all series are z-normalized in this case, since
#'       the resulting centroids will also have this normalization. See
#'       \code{\link{shape_extraction}} for more details.
#'     \item "dba": DTW Barycenter Averaging. See \code{\link{DBA}} for more details.
#'     \item "pam": Partition around medoids (PAM). This basically means that the cluster centroids
#'       are always one of the time series in the data. In this case, the distance matrix can be
#'       pre-computed once using all time series in the data and then re-used at each iteration. It
#'       usually saves overhead overall.
#'     \item "fcm": Fuzzy c-means. Only supported for fuzzy clustering and used by default in that
#'       case.
#'     \item "fcmdd": Fuzzy c-medoids. Only supported for fuzzy clustering. It \strong{always}
#'       precomputes the whole cross-distance matrix.
#'   }
#'
#'   These check for the special cases where parallelization might be desired. Note that only
#'   \code{shape}, \code{dba} and \code{pam} support series of different length. Also note that, for
#'   \code{shape} and \code{dba}, this support has a caveat: the final centroids' length will depend
#'   on the length of those series that were randomly chosen at the beginning of the clustering
#'   algorithm. For example, if the series in the dataset have a length of either 10 or 15, 2
#'   clusters are desired, and the initial choice selects two series with length of 10, the final
#'   centroids will have this same length.
#'
#'   As special cases, if hierarchical or tadpole clustering is used, you can provide a centroid
#'   function that takes a list of series as only input and returns a single centroid series. These
#'   centroids are returned in the \code{centroids} slot. By default, a type of PAM centroid
#'   function is used.
#'
#' @section Distance Measures:
#'
#'   The distance measure to be used with partitional, hierarchical and fuzzy clustering can be
#'   modified with the \code{distance} parameter. The supported option is to provide a string, which
#'   must represent a compatible distance registered with \code{proxy}'s \code{\link[proxy]{dist}}.
#'   Registration is done via \code{\link[proxy]{pr_DB}}, and extra parameters can be provided in
#'   \code{args$dist}.
#'
#'   Note that you are free to create your own distance functions and register them. Optionally, you
#'   can use one of the following custom implementations (all registered with \code{proxy}):
#'
#'   \itemize{
#'     \item \code{"dtw"}: DTW, optionally with a Sakoe-Chiba/Slanted-band constraint*.
#'     \item \code{"dtw2"}: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint*.
#'       Read details below.
#'     \item \code{"dtw_basic"}: A custom version of DTW with less functionality, but slightly
#'       faster. See \code{\link{dtw_basic}}.
#'     \item \code{"dtw_lb"}: DTW with L1 or L2 norm* and optionally a Sakoe-Chiba constraint*. Some
#'       computations are avoided by first estimating the distance matrix with Lemire's lower bound
#'       and then iteratively refining with DTW. See \code{\link{dtw_lb}}. Not suitable for
#'       \code{pam.precompute}* = \code{TRUE}.
#'     \item \code{"lbk"}: Keogh's lower bound with either L1 or L2 norm* for the Sakoe-Chiba
#'       constraint*.
#'     \item \code{"lbi"}: Lemire's lower bound with either L1 or L2 norm* for the Sakoe-Chiba
#'       constraint*.
#'     \item \code{"sbd"}: Shape-based distance. See \code{\link{SBD}} for more details.
#'     \item \code{"gak"}: Global alignment kernels. See \code{\link{GAK}} for more details.
#'   }
#'
#'   DTW2 is done with \code{\link[dtw]{dtw}}, but it differs from the result you would obtain if
#'   you specify \code{L2} as \code{dist.method}: with \code{DTW2}, pointwise distances (the local
#'   cost matrix) are calculated with \code{L1} norm, \emph{each} element of the matrix is squared
#'   and the result is fed into \code{\link[dtw]{dtw}}, which finds the optimum warping path. The
#'   square root of the resulting distance is \emph{then} computed. See \code{\link{dtw2}}.
#'
#'   Only \code{dtw}, \code{dtw2}, \code{sbd} and \code{gak} support series of different length. The
#'   lower bounds are probably unsuitable for direct clustering unless series are very easily
#'   distinguishable.
#'
#'   If you know that the distance function is symmetric, and you use a hierarchical algorithm, or a
#'   partitional algorithm with PAM centroids and \code{pam.precompute} = \code{TRUE}, some time
#'   can be saved by calculating only half the distance matrix. Therefore, consider setting the
#'   symmetric control parameter to \code{TRUE} if this is the case.
#'
#' @section Preprocessing:
#'
#'   It is strongly advised to use z-normalization in case of \code{centroid = "shape"}, because the
#'   resulting series have this normalization (see \code{\link{shape_extraction}}). Therefore,
#'   \code{\link{zscore}} is the default in this case. The user can, however, specify a custom
#'   function that performs any transformation on the data, but the user must make sure that the
#'   format stays consistent, i.e. a list of time series.
#'
#'   Setting to \code{NULL} means no preprocessing (except for \code{centroid = "shape"}). A
#'   provided function will receive the data as first argument, followed by the contents of
#'   \code{args$preproc} that match its formal arguments.
#'
#'   It is convenient to provide this function if you're planning on using the
#'   \code{\link[stats]{predict}} generic.
#'
#' @section Repetitions:
#'
#'   Due to their stochastic nature, partitional clustering is usually repeated* several times with
#'   different random seeds to allow for different starting points. This function uses
#'   \code{\link[rngtools]{RNGseq}} to obtain different seed streams for each repetition, utilizing
#'   the \code{seed} parameter (if provided) to initialize it. If more than one repetition is made,
#'   the streams are returned in an attribute called \code{rng}.
#'
#'   Multiple values of \code{k} can also be provided to get different partitions using any
#'   \code{type} of clustering.
#'
#'   Repetitions are greatly optimized when PAM centroids are used and the whole distance matrix is
#'   precomputed, since said matrix is reused for every repetition, and can be comptued in parallel
#'   (see Parallel section).
#'
#' @section Parallel Computing:
#'
#'   Please note that running tasks in parallel does \strong{not} guarantee faster computations. The
#'   overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#'   The user can register a parallel backend, for eample with the \code{doParallel} package, in
#'   order to do the repetitions in parallel, as well as distance and some centroid calculations.
#'
#'   Unless each repetition requires a few seconds, parallel computing probably isn't worth it. As
#'   such, I would only use this feature with \code{shape} and \code{DBA} centroids, or an expensive
#'   distance function like \code{DTW}.
#'
#'   If you register a parallel backend, the function will also try to do the calculation of the
#'   distance matrices in parallel. This should work with any function registered with
#'   \code{\link[proxy]{dist}} via \code{\link[proxy]{pr_DB}} whose \code{loop} flag is set to
#'   \code{TRUE}. If the function requires special packages to be loaded, provide their names in the
#'   \code{packages}* slot of \code{control}. Note that "dtwclust" is always loaded in each parallel
#'   worker, so that doesn't need to be included. Alternatively, you may want to pre-load
#'   \code{dtwclust} in each worker with \code{\link[parallel]{clusterEvalQ}}.
#'
#'   In case of multiple repetitions, each worker gets a repetition task. Otherwise, the tasks
#'   (which can be a distance matrix or a centroid calculation) are usually divided into chunks
#'   according to the number of workers available.
#'
#' @section Notes:
#'
#'   The lower bounds are defined only for time series of equal length. \code{DTW} and \code{DTW2}
#'   don't require this, but they are much slower to compute.
#'
#'   The lower bounds are \strong{not} symmetric, and \code{DTW} is not symmetric in general.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references.
#'
#' @seealso
#'
#' \code{\link{TSClusters-class}}, \code{\link{tsclusters-methods}},
#' \code{\link{tsclustFamily-class}}.
#'
#' @example inst/tsclust-examples.R
#'
tsclust <- function(series = NULL, type = "partitional", k = 2L, ...,
                    preproc = NULL, distance = "dtw_basic",
                    centroid = ifelse(type == "fuzzy", "fcm", "pam"),
                    control = do.call(paste0(type, "_control"), list()),
                    args = tsclust_args(),
                    seed = NULL, trace = FALSE)
{
    ## =============================================================================================
    ## Start
    ## =============================================================================================

    tic <- proc.time()

    set.seed(seed)

    if (is.null(series))
        stop("No data provided")

    type <- match.arg(type, c("partitional", "hierarchical", "tadpole", "fuzzy"))

    ## coerce to list if necessary
    series <- any2list(series)

    if (any(k < 2L))
        stop("At least two clusters must be defined")
    if (any(k > length(series)))
        stop("Cannot have more clusters than series in the data")
    if (!is.list(control))
        stop("Invalid control argument")

    MYCALL <- match.call(expand.dots = TRUE)

    dots <- list(...)

    ## ---------------------------------------------------------------------------------------------
    ## Preprocess
    ## ---------------------------------------------------------------------------------------------

    if (!is.null(preproc) && is.function(preproc)) {
        series <- do.call(preproc,
                          enlist(series,
                                 dots = subset_dots(args$preproc, preproc)))

        preproc_char <- as.character(substitute(preproc))[1L]

    } else if (type == "partitional" && is.character(centroid) && centroid == "shape") {
        preproc <- zscore
        preproc_char <- "zscore"
        series <- do.call(zscore,
                          enlist(series,
                                 dots = subset_dots(args$preproc, zscore)))

    } else if (is.null(preproc)) {
        preproc <- function(x, ...) { x }
        preproc_char <- "none"

    } else stop("Invalid preprocessing")

    check_consistency(series, "vltslist")

    ## ---------------------------------------------------------------------------------------------
    ## Further options
    ## ---------------------------------------------------------------------------------------------

    diff_lengths <- different_lengths(series)

    check_consistency(distance, "dist", trace = trace, Lengths = diff_lengths, silent = FALSE)

    if (type %in% c("partitional", "hierarchical")) {
        ## symmetric versions of dtw that I know of
        ## unconstrained and with symmetric1/symmetric2 is always symmetric, regardless of lengths
        ## constrained and same lengths with symmetric1/symmetric2 is also symmetric
        symmetric_pattern <- is.null(dots$step.pattern) ||
            identical(dots$step.pattern, symmetric1) ||
            identical(dots$step.pattern, symmetric2)

        if (tolower(distance) %in% c("dtw", "dtw2", "dtw_basic")) {
            control$symmetric <- symmetric_pattern && (is.null(dots$window.size) || !diff_lengths)

        } else if (tolower(distance) %in% c("lbk", "lbi")) {
            control$symmetric <- FALSE

        } else if (tolower(distance) %in% c("sbd", "gak")) {
            control$symmetric <- TRUE
        }
    }

    RET <-
        switch(type,
               partitional =, fuzzy = {
                   ## ==============================================================================
                   ## Partitional or fuzzy
                   ## ==============================================================================

                   if (!inherits(control, "PtCtrl") && !inherits(control, "FzCtrl"))
                       stop("Invalid control provided")

                   nrep <- if (is.null(control$nrep)) 1L else control$nrep

                   if (is.character(centroid)) {
                       if (type == "fuzzy")
                           centroid <- match.arg(centroid, c("fcm", "fcmdd"))
                       else
                           centroid <- match.arg(centroid, c("mean", "median", "shape", "dba", "pam"))

                       ## replace any given distmat if centroid not "pam" or "fcmdd"
                       if (!(centroid %in% c("pam", "fcmdd")))
                           distmat <- NULL
                       else
                           distmat <- control$distmat

                   } else {
                       distmat <- NULL
                   }

                   if (diff_lengths) {
                       if (type == "fuzzy" && is.character(centroid) && centroid == "fcm")
                           stop("Fuzzy c-means does not support series with different length.")

                       check_consistency(centroid, "cent", trace = trace)
                   }

                   ## ------------------------------------------------------------------------------
                   ## Family creation, see initialization in tsclusters-methods.R
                   ## ------------------------------------------------------------------------------

                   family <- new("tsclustFamily",
                                 dist = distance,
                                 allcent = centroid,
                                 preproc = preproc,
                                 distmat = distmat,
                                 control = control,
                                 fuzzy = isTRUE(type == "fuzzy"))

                   if (!all(c("x", "cl_id", "k", "cent", "cl_old") %in% names(formals(family@allcent))))
                       stop("The provided centroid function must have at least the following ",
                            "arguments with the shown names:\n\t",
                            paste(c("x", "cl_id", "k", "cent", "cl_old"), collapse = ", "))

                   cent_char <- as.character(substitute(centroid))[1L]

                   ## ------------------------------------------------------------------------------
                   ## PAM precompute?
                   ## ------------------------------------------------------------------------------

                   ## for a check near the end, changed if appropriate
                   distmat_provided <- FALSE

                   ## precompute distance matrix?
                   if (cent_char %in% c("pam", "fcmdd")) {
                       ## check if distmat was not provided and should be precomputed
                       if (!is.null(distmat)) {
                           if ( nrow(distmat) != length(series) || ncol(distmat) != length(series) )
                               stop("Dimensions of provided cross-distance matrix don't correspond ",
                                    "to length of provided data")

                           ## distmat was provided in call
                           distmat_provided <- TRUE

                           if (trace) cat("\n\tDistance matrix provided...\n\n")

                       } else if (isTRUE(control$pam.precompute) || cent_char == "fcmdd") {
                           if (tolower(distance) == "dtw_lb")
                               warning("Using dtw_lb with control$pam.precompute = TRUE is not ",
                                       "advised.")

                           if (trace) cat("\n\tPrecomputing distance matrix...\n\n")

                           distmat <- do.call(family@dist,
                                              enlist(x = series,
                                                     centroids = NULL,
                                                     dots = args$dist))

                           ## Redefine new distmat
                           environment(family@dist)$control$distmat <- distmat
                           environment(family@allcent)$distmat <- distmat

                           gc(FALSE)
                       }
                   }

                   ## ------------------------------------------------------------------------------
                   ## Cluster
                   ## ------------------------------------------------------------------------------

                   if (length(k) == 1L && nrep == 1L) {
                       rngtools::setRNG(rngtools::RNGseq(1L, seed = seed, simplify = TRUE))

                       ## Just one repetition
                       pc.list <- list(do.call(pfclust,
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

                       rng <- rngtools::RNGseq(length(k) * nrep,
                                               seed = seed, simplify = FALSE)
                       rng0 <- lapply(parallel::splitIndices(length(rng), length(k)),
                                      function(i) rng[i])

                       ## if %do% is used, the outer loop replaces value of k in this envir
                       k0 <- k
                       comb0 <- if (nrep > 1L) c else list

                       i <- integer() # CHECK complains about non-initialization now

                       pc.list <- foreach(k = k0, rng = rng0, .combine = comb0, .multicombine = TRUE,
                                          .packages = control$packages, .export = export) %:%
                           foreach(i = 1L:nrep, .combine = list, .multicombine = TRUE,
                                   .packages = control$packages, .export = export) %op%
                                   {
                                       if (trace) message("Repetition ", i, " for k = ", k)

                                       rngtools::setRNG(rng[[i]])

                                       if (!check_consistency(dist_entry$names[1], "dist"))
                                           do.call(proxy::pr_DB$set_entry, dist_entry)

                                       pc <- do.call(pfclust,
                                                     enlist(x = series,
                                                            k = k,
                                                            family = family,
                                                            control = control,
                                                            fuzzy = isTRUE(type == "fuzzy"),
                                                            cent = cent_char,
                                                            trace = trace,
                                                            args = args))

                                       gc(FALSE)

                                       pc
                                   }
                   }

                   ## ------------------------------------------------------------------------------
                   ## Prepare results
                   ## ------------------------------------------------------------------------------

                   ## Replace distmat with NULL so that, if the distance function is called again,
                   ## it won't subset it
                   eval(quote(control$distmat <- NULL), environment(family@dist))
                   eval(quote(distmat <- NULL), environment(family@allcent))

                   ## If distmat was provided, let it be shown in the results
                   if (distmat_provided) {
                       if (is.null(attr(distmat, "method")))
                           distance <- "unknown"
                       else
                           distance <- attr(distmat, "method")
                   }

                   ## Create objects
                   RET <- lapply(pc.list, function(pc) {
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

                   if (!inherits(RET, "TSClusters") && length(RET) == 1L) RET <- RET[[1L]]

                   RET
               },

               hierarchical = {
                   ## ==============================================================================
                   ## Hierarchical
                   ## ==============================================================================

                   if (!inherits(control, "HcCtrl")) stop("Invalid control provided")

                   method <- control$method
                   distmat <- control$distmat

                   if (tolower(distance) == "dtw_lb")
                       warning("Using dtw_lb with hierarchical clustering is not advised.")

                   ## ------------------------------------------------------------------------------
                   ## Calculate distance matrix
                   ## ------------------------------------------------------------------------------

                   ## Take advantage of the function I defined for the partitional methods
                   ## Which can do calculations in parallel if appropriate
                   distfun <- ddist2(distance = distance,
                                     control = control)

                   if (!is.null(distmat)) {
                       if (nrow(distmat) != length(series) || ncol(distmat) != length(series))
                           stop("Dimensions of provided cross-distance matrix don't correspond to ",
                                "length of provided data")

                       if (trace)
                           cat("\n\tDistance matrix provided...\n")

                       if (is.null(attr(distmat, "method")))
                           distance <- "unknown"
                       else
                           distance <- attr(distmat, "method")

                   } else {
                       if (trace)
                           cat("\n\tCalculating distance matrix...\n")

                       ## single argument is to calculate whole distance matrix
                       distmat <- do.call(distfun, enlist(x = series,
                                                          centroids = NULL,
                                                          dots = args$dist))
                   }

                   ## ------------------------------------------------------------------------------
                   ## Cluster
                   ## ------------------------------------------------------------------------------

                   if (trace)
                       cat("\n\tPerforming hierarchical clustering...\n\n")

                   if (is.character(method)) {
                       ## Using hclust
                       hc <- lapply(method, function(method) {
                           stats::hclust(stats::as.dist(distmat), method, members = dots$members)
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

                       method <- attr(method, "name")
                   }

                   ## Invalid centroid specifier provided?
                   if (!missing(centroid) && !is.function(centroid))
                       warning("The 'centroid' argument was provided but it wasn't a function, ",
                               "so it was ignored.")
                   if (!is.function(centroid))
                       centroid <- NA

                   ## ------------------------------------------------------------------------------
                   ## Prepare results
                   ## ------------------------------------------------------------------------------

                   RET <- lapply(k, function(k) {
                       lapply(hc, function(hc) {
                           ## cutree and corresponding centroids
                           cluster <- stats::cutree(stats::as.hclust(hc), k)

                           if (is.function(centroid)) {
                               allcent <- centroid
                               centroids <- lapply(1L:k, function(kcent) {
                                   do.call(centroid,
                                           enlist(series[cluster == kcent],
                                                  dots = subset_dots(args$cent, centroid)))
                               })
                               cent_char <- as.character(substitute(centroid))[1L]

                           } else {
                               allcent <- function() {}

                               centroids <- sapply(1L:k, function(kcent) {
                                   id_k <- cluster == kcent

                                   d_sub <- distmat[id_k, id_k, drop = FALSE]

                                   id_centroid <- which.min(apply(d_sub, 1L, sum))

                                   which(id_k)[id_centroid]
                               })

                               centroids <- series[centroids]
                               cent_char <- "PAM (Hierarchical)"
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

                   if (!inherits(RET, "TSClusters") && length(RET) == 1L) RET <- RET[[1L]]

                   RET
               },

               tadpole = {
                   ## ==============================================================================
                   ## TADPole
                   ## ==============================================================================

                   if (!inherits(control, "TpCtrl")) stop("Invalid control provided")
                   if (!missing(distance)) warning("The distance argument is ignored for TADPole.")

                   ## ------------------------------------------------------------------------------
                   ## Cluster
                   ## ------------------------------------------------------------------------------

                   ## mainly for predict generic
                   distfun <- ddist2("dtw_lb", control = control)

                   ## Invalid centroid specifier provided?
                   if (!missing(centroid) && !is.function(centroid))
                       warning("The 'centroid' argument was provided but it wasn't a function, ",
                               "so it was ignored.")

                   if (is.function(centroid))
                       cent_char <- as.character(substitute(centroid))
                   else
                       cent_char <- "PAM (TADPole)"

                   ## for family@dist
                   args$dist$window.size <- control$window.size
                   args$dist$norm <- "L2"

                   ## seeds
                   rng <- rngtools::RNGseq(length(k), seed = seed, simplify = FALSE)

                   RET <- foreach(k = k, rng = rng,
                                  .combine = list, .multicombine = TRUE,
                                  .packages = "dtwclust", .export = "enlist") %op%
                                  {
                                      rngtools::setRNG(rng)

                                      if (trace) cat("\nEntering TADPole...\n\n")

                                      R <- TADPole(series, k = k,
                                                   dc = control$dc,
                                                   window.size = control$window.size,
                                                   lb = control$lb)

                                      if (trace) {
                                          cat("TADPole completed, pruning percentage = ",
                                              formatC(100 - R$distCalcPercentage,
                                                      digits = 3L,
                                                      width = -1L,
                                                      format = "fg"),
                                              "%\n\n",
                                              sep = "")
                                      }

                                      ## -----------------------------------------------------------
                                      ## Prepare results
                                      ## -----------------------------------------------------------

                                      if (is.function(centroid)) {
                                          allcent <- centroid
                                          centroids <- lapply(1L:k, function(kcent) {
                                              centroid(series[R$cl == kcent])
                                          })

                                      } else {
                                          allcent <- function(x) {}
                                          centroids <- series[R$centroids]
                                      }

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
                                  }
               })

    ## =============================================================================================
    ## Finish
    ## =============================================================================================

    toc <- proc.time() - tic

    if (inherits(RET, "TSClusters")) {
        RET@proctime <- toc

    } else {
        RET <- lapply(RET, function(ret) {
            ret@proctime <- toc

            ret
        })
    }

    if (type %in% c("partitional", "fuzzy") && (nrep > 1L || length(k) > 1L))
        attr(RET, "rng") <- unlist(rng0, recursive = FALSE, use.names = FALSE)

    if (trace)
        cat("\tElapsed time is", toc["elapsed"], "seconds.\n\n")

    RET
}
