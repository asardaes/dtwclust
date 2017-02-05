#' Time series clustering
#'
#' This is the new experimental main function to perform time series clustering. It should provide
#' the same functionality as \code{\link{dtwclust}}, but it is hopefully more coherent in general.
#' \strong{For now, it is subject to change}. Feedback is appreciated at
#' \url{https://github.com/asardaes/dtwclust/issues}.
#'
#' @export
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
                           assign("distmat", distmat, environment(family@dist))
                           assign("distmat", distmat, environment(family@allcent))

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

                   if (class(RET) != "dtwclust" && length(RET) == 1L) RET <- RET[[1L]]

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
                                     control = control,
                                     distmat = NULL)

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
                   distfun <- ddist2("dtw_lb", control = control, distmat = NULL)

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
