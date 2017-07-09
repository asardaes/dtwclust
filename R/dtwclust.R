#' Time series clustering
#'
#' This is the **deprecated** function to perform time series clustering. See [tsclust()] for the
#' new interface.
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
#' `DTW2` and `GAK` support such series, which means only partitional and hierarchical procedures
#' using those distances will work. You can of course create your own custom distances. All included
#' centroid functions should work with the aforementioned format, although `shape` is **not**
#' recommended. Note that the `plot` method will simply append all dimensions (columns) one after
#' the other.
#'
#' Several parameters can be adjusted with the `control` argument. See [dtwclustControl]. In the
#' following sections, elements marked with an asterisk (*) are those that can be adjusted with this
#' argument.
#'
#' @return
#'
#' An object with formal class [dtwclust-class].
#'
#' If `control@nrep > 1` and a partitional procedure is used, `length(method)` `> 1` and
#' hierarchical procedures are used, or `length(k)` `>` `1`, a list of objects is returned.
#'
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references.
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
    message("This FUNCTION is now deprecated and will be eventually removed.\n",
            "Please use 'dtwclust::tsclust' instead. See ?'dtwclust-deprecated'")

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

    if (is.null(control) || is.list(control)) control <- methods::as(control, "dtwclustControl")
    else if (class(control) != "dtwclustControl") stop("Invalid control argument")
    else methods::validObject(control)

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

    check_consistency(distance, "dist", trace = control@trace, diff_lengths = diff_lengths, silent = FALSE)

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

            if (is.character(centroid) && (centroid %in% centroids_included) && !(centroid %in% centroids_difflength))
                stop("Only the following centroids are supported for series with different lengths:\n\t",
                     paste(centroids_difflength, collapse = "\t"))
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
