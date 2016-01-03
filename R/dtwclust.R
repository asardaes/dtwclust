#' Time series clustering
#'
#' This function uses the DTW distance and related techniques to cluster time series. It also supports custom
#' distances, especially through \code{\link[proxy]{dist}}, as well as custom centroid functions.
#' For now, all series must be univariate.
#'
#' Partitional algorithms use custom functions inspired by \code{\link[flexclust]{kcca}}. Hierarchical algorithms use
#' the \code{\link[stats]{hclust}} function. The \code{tadpole} algorithm uses the \code{\link{TADPole}} function.
#'
#' The \code{data} may be a matrix or a list, but the matrix will be coerced to a list row-wise. A matrix requires
#' that all time series have equal lengths. If the lengths vary slightly between time series, reinterpolating them to
#' a common length is most likely an acceptable approach (Ratanamahatana and Keogh, 2004). See the examples.
#'
#' Several parameters can be adjusted with the \code{control} argument. See \code{\link{dtwclustControl}}.
#'
#' @section Distance:
#'
#' The supported option is to provide a string. The string can represent a compatible registered distance of
#' \code{proxy}'s \code{\link[proxy]{dist}}. Extra parameters can be provided in \code{...}.
#'
#' It can be a string of one of the following custom implementations (all registered with \code{proxy}):
#'
#' \itemize{
#'   \item \code{"dtw"}: DTW with L1 norm and optionally a Sakoe-Chiba/Slanted-band constraint.
#'   \item \code{"dtw2"}: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint.
#'   \item \code{"dtw_lb"}: DTW with L1 or L2 norm and optionally a Sakoe-Chiba constraint. Some computations
#'   are avoided by first estimating the distance matrix with Lemire's lower bound and then iteratively
#'   refining with DTW. See \code{\link{dtw_lb}}. Not suitable for \code{pam.precompute} = \code{TRUE}
#'   (see \code{\link{dtwclustControl}}).
#'   \item \code{"lbk"}: Keogh's lower bound with either L1 or L2 norm for the Sakoe-Chiba constraint.
#'   \item \code{"lbi"}: Lemire's lower bound with either L1 or L2 norm for the Sakoe-Chiba constraint.
#'   \item \code{"sbd"}: Shape-based distance. Each series is z-normalized in this case. As a result,
#'   the cluster centers (for partitional methods) are also z-normalized. See \code{\link{SBD}} for more
#'   details.
#' }
#'
#' Note that only \code{dtw}, \code{dtw2} and \code{sbd} support series of different lengths. The lower bounds
#' are probably unsuitable for direct clustering unless series are very easily distinguishable.
#'
#' If you create your own distance and also register it with \code{proxy} (see \code{\link[proxy]{pr_DB}}),
#' and it includes the ellipsis (\code{...}) in its definition, it will receive
#' the following parameters (see \code{\link{dtwclustControl}}):
#'
#' \itemize{
#'   \item \code{window.type}: Either \code{"none"} for a \code{NULL} \code{window.size}, or \code{"slantedband"}
#'   otherwise
#'   \item \code{window.size}: The provided window size
#'   \item \code{norm}: The provided desired norm
#'   \item \code{...}: Any additional parameters provided in the original call's ellipsis
#' }
#'
#' Whether the function makes use of them or not, is up to you. See the examples.
#'
#' If you know that the distance function is symmetric, consider setting the corresponding slot in
#' \code{\link{dtwclustControl}} to \code{TRUE} to save some time in case of \code{pam.precompute}
#' = \code{TRUE} or hierarchical procedures.
#'
#' @section Centroid:
#'
#' In the case of partitional algorithms, a suitable function should calculate the cluster centers at
#' every iteration. In this case, the centers are themselves time series.
#'
#' If a custom function is provided, it will receive the following inputs in the shown order (examples
#' shown in parenthesis):
#'
#' \itemize{
#'   \item The \emph{whole} data list (\code{list(ts1, ts2, ts3)})
#'   \item A numeric vector with length equal to the number of series in \code{data}, indicating which
#'   cluster a series belongs to (\code{c(1L, 2L, 2L)})
#'   \item The desired number of total clusters (\code{2L})
#'   \item The current centers in order, in a list (\code{list(center1, center2)})
#'   \item The membership vector of the \emph{previous} iteration (\code{c(1L, 1L, 2L)})
#'   \item The elements of \code{...}
#' }
#'
#' Therefore, the function should always include the ellipsis \code{...} in its definition.
#'
#' The other option is to provide a character string for the custom implementations. The following options
#' are available:
#'
#' \itemize{
#'   \item "mean": The average along each dimension. In other words, the average of all \eqn{x^j_i}
#'   among the \eqn{j} series that belong to the same cluster for all time points \eqn{t_i}.
#'   \item "median": The median along each dimension. Similar to mean.
#'   \item "shape": Shape averaging. See \code{\link{shape_extraction}} for more details.
#'   \item "dba": DTW Barycenter Averaging. See \code{\link{DBA}} for more details.
#'   \item: "pam": Partition around medoids. This basically means that the cluster centers are always
#'   one of the time series in the data. In this case, the distance matrix can be pre-computed once using all
#'   time series in the data and then re-used at each iteration. It usually saves overhead overall.
#' }
#'
#' These check for the special cases where parallelization is desired.
#'
#' If a cluster becomes empty, the functions reinitialize a new cluster randomly. However, if you see
#' that the algorithm doesn't converge or the overall distance sum increases, this probably means that
#' the chosen value of \code{k} is too large, or the chosen distance function is not able to assess
#' similarity effectively. The random reinitialization attempts to enforce a certain number of clusters,
#' but can result in instability in the aforementioned cases.
#'
#' Note that only \code{shape}, \code{dba} and \code{pam} support series of different lengths.
#'
#' @section Sakoe-Chiba Constraint:
#'
#' A global constraint to speed up the DTW calculation is the Sakoe-Chiba band (Sakoe and Chiba, 1978). To
#' use it, a window width must be defined via the \code{window.size} slot in \code{control} (see
#' \code{\link{dtwclustControl}}).
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10*2 + 1 = 21} observations falling within
#' the window.
#'
#' The computations actually use a \code{slantedband} window, which is equivalent to the Sakoe-Chiba one
#' if series have equal lengths, and stays along the diagonal of the local cost matrix if series have
#' different lengths.
#'
#' @section Preprocessing:
#'
#' It is strongly advised to use z-normalization in case of \code{centroid = "shape"}, because the resulting
#' series have this normalization (see \code{\link{shape_extraction}}). Therefore, \code{\link{zscore}} is the
#' default in this case. The user can, however, specify a
#' custom function that performs any transformation on the data, but the user must make sure that the format
#' stays consistent, i.e. a list of time series.
#'
#' Setting to \code{NULL} means no preprocessing (except for \code{centroid = "shape"}).
#'
#' A provided function will receive the data as first and only argument.
#'
#' It is convenient to provide this function if you're planning on using the \code{\link[stats]{predict}} generic.
#'
#' @section Repetitions:
#'
#' Due to their stochastic nature, partitional clustering is usually repeated several times with different random
#' seeds to allow for different starting points. This function can now run several repetitions by using the
#' \code{doRNG} package. This package ensures that each repetition uses a statistically independent random sequence.
#'
#' If more than one repetition is made, the \code{\link[doRNG]{\%dorng\%}} operator is used. If provided, the
#' \code{seed} parameter is used to initialize it. The different seed sequences used are returned in the
#' \code{rng} attribute in such cases.
#'
#' Repetitions are greatly optimized when PAM centroids are used and the whole distance matrix is precomputed,
#' since said matrix is reused for every repetition, and can be comptued in parallel (see next section).
#'
#' The number of desired repetitions is provided through the \code{control} parameter. See \code{\link{dtwclustControl}}.
#'
#' @section Parallel Computing:
#'
#' Please note that running tasks in parallel does \strong{not} guarantee faster computations.
#' The overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#' The user can register a parallel backend with the \code{doParallel} package (and possibly any other
#' package compatible with \code{foreach}'s \code{\%dopar\%} operator) in order to do the
#' repetitions in parallel, as well as distance and some centroid calculations (see the examples).
#' \code{\link{TADPole}} and \code{\link{DBA}} also take advantage of parallel support.
#'
#' Unless each repetition requires a few seconds, parallel computing probably isn't worth it. As such, I would only
#' use this feature with \code{shape} and \code{DBA} centroids, or an expensive distance function like \code{DTW}
#' (see the next paragraph).
#'
#' If you register a parallel backend, the function will also try to do the calculation of the distance
#' matrices in parallel. This should work with any function registered with \code{\link[proxy]{dist}} via
#' \code{\link[proxy]{pr_DB}} whose \code{loop} flag is set to \code{TRUE}. If the function requires special packages
#' to be loaded, provide their names in the \code{packages} slot of \code{control} (see \code{\link{dtwclustControl}}).
#' In addition, "dtwclust" is always loaded in each parallel worker, so that doesn't need to be included.
#' Alternatively, you may want to pre-load \code{dtwclust} in each worker with \code{\link[parallel]{clusterEvalQ}}.
#'
#' Note that, by default, if a parallel backend is registered, multiple repetitions are to be performed (partitional clustering,
#' \code{control@nrep} \code{>=} 1)
#' AND \code{centroid} \code{!=} \code{"pam"}, each parallel worker will get a repetition task, but any calculations
#' within each worker will be done sequentially. Load balance for such a scenario should be fine as long as the
#' \code{control@nrep} parameter is greater than or equal to the number of parallel workers.
#' If you believe your task would benefit more from parallelization within each repetition,
#' consider registering the parallel backend and calling \code{dtwclust} several times sequentially,
#' with \code{control@nrep = 1} and different \code{seeds}.
#'
#' @section Notes:
#'
#' The lower bounds are defined only for time series of equal lengths. \code{DTW} and \code{DTW2}
#' don't require this, but they are much slower to compute.
#'
#' The lower bounds are \strong{not} symmetrical, and \code{DTW} is only symmetrical if series are of equal
#' lengths.
#'
#' Specifying \code{distance = "sbd"} and \code{centroid = "shape"} is equivalent to the k-Shape algorithm
#' (Papparizos and Gravano, 2015). See \code{\link{SBD}} and \code{\link{shape_extraction}} for more info.
#'
#' Hierarchical (conditional on \code{method}) and TADpole (conditional on \code{dc}) procedures are deterministic.
#'
#' DTW2 is done with \code{\link[dtw]{dtw}},but it differs from the result you would obtain if you specify
#' \code{L2} as \code{dist.method}: with \code{DTW2}, pointwise distances (the local cost matrix) are calculated with
#' \code{L1} norm, \emph{each} element of the matrix is squared and the result is fed into
#' \code{\link[dtw]{dtw}}, which finds the optimum warping path. The square root of the resulting
#' distance is \emph{then} computed.
#'
#' @references
#'
#' Sakoe H and Chiba S (1978). ``Dynamic programming algorithm optimization for spoken word recognition.'' \emph{Acoustics, Speech
#' and
#' Signal Processing, IEEE Transactions on}, \strong{26}(1), pp. 43-49. ISSN 0096-3518, \url{http://doi.org/10.1109/TASSP.1978.1163055}.
#'
#' Ratanamahatana A and Keogh E (2004). ``Everything you know about dynamic time warping is wrong.'' In \emph{3rd Workshop on Mining
#' Temporal and Sequential Data, in conjunction with 10th ACM SIGKDD Int. Conf. Knowledge Discovery and Data Mining (KDD-2004),
#' Seattle, WA}.
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.'' In \emph{Proceedings of the 2015
#' ACM SIGMOD International Conference on Management of Data}, series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{
#' http://doi.org/10.1145/2723372.2737793}.
#'
#' @example inst/dtwclust-examples.R
#'
#' @seealso
#'
#' Please check the brief description in \code{\link{dtwclust-package}}.
#'
#' Type \code{news(package = "dtwclust")} to see what changed.
#'
#' Additionally: \code{\link{dtwclust-methods}}, \code{\link{dtwclust-class}}, \code{\link{dtwclustControl}},
#' \code{\link{dtwclustFamily}}.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @param data A list where each element is a time series, or a numeric matrix (it will be coerced to a list
#' row-wise).
#' @param type What type of clustering method to use: \code{partitional}, \code{hierarchical} or \code{tadpole}.
#' @param k Numer of desired clusters in partitional methods. For hierarchical methods, the
#' \code{\link[stats]{cutree}} function is called with this value of \code{k} and the result is returned in the
#' \code{cluster} slot of the \code{dtwclust} object.
#' @param method One or more linkage methods to use in hierarchical procedures. See \code{\link[stats]{hclust}}.
#' You can provide a character vector to compute different hierarchical cluster structures in one go, or
#' specify \code{method} = "all" to use all the available ones.
#' @param distance A supported distance from \code{\link[proxy]{pr_DB}} (see Distance section). Ignored for
#' \code{type} = \code{"tadpole"}.
#' @param centroid Either a supported string or an appropriate function to calculate centroids
#' when using partitional methods (see Centroid section).
#' @param preproc Function to preprocess data. Defaults to \code{zscore} \emph{only} if \code{centroid}
#' \code{=} \code{"shape"}, but will be replaced by a custom function if provided. See Preprocessing section.
#' @param dc Cutoff distance for the \code{\link{TADPole}} algorithm.
#' @param control NEW: Named list of parameters for clustering algorithms. See \code{\link{dtwclustControl}}.
#' \code{NULL} means defaults.
#' @param seed Random seed for reproducibility of partitional algorithms.
#' @param distmat If a cross-distance matrix is already available, it can be provided here so it's re-used. Only relevant if
#' \code{centroid} = "pam" or \code{type} = "hierarchical". See examples.
#' @param ... Additional arguments to pass to \code{\link[proxy]{dist}} or a custom function. For now, old parameters
#' can also be provided here. See \code{\link{dtwclustControl}}.
#'
#' @return An object with formal class \code{\link{dtwclust-class}}.
#' If \code{control@nrep > 1} and a partitional procedure is used, or \code{length(method)} \code{> 1} and hierarchical
#' procedures are used, a list of objects is returned.
#'
#' @export
#'

dtwclust <- function(data = NULL, type = "partitional", k = 2L, method = "average",
                     distance = "dtw", centroid = "pam", preproc = NULL,
                     dc = NULL, control = NULL, seed = NULL, distmat = NULL,
                     ...)
{
     ## =================================================================================================================
     ## Start
     ## =================================================================================================================

     tic <- proc.time()

     if (is.null(data))
          stop("No data provided")

     type <- match.arg(type, c("partitional", "hierarchical", "tadpole"))

     ## coerce to list if necessary
     data <- consistency_check(data, "tsmat")

     MYCALL <- match.call(expand.dots = TRUE)

     ## ----------------------------------------------------------------------------------------------------------
     ## Control parameters
     ## ----------------------------------------------------------------------------------------------------------

     if (is.null(control) || is.list(control)) {
          control <- as(control, "dtwclustControl")

     } else if (class(control) != "dtwclustControl" || !validObject(control)) {
          stop("Invalid control argument") # validObject should generate an error anyway
     }

     ## ----------------------------------------------------------------------------------------------------------
     ## Backwards compatibility, check for old formal arguments in '...'
     ## I might leave this here to prevent duplicate matching anyway
     ## ----------------------------------------------------------------------------------------------------------

     dots <- list(...)

     args <- names(MYCALL[-1])
     id_oldargs <- which(args %in% c(slotNames("dtwclustControl"), "reps"))
     if (length(id_oldargs) > 0) {
          message("The 'control' argument replaced many parameters that were dropped from this function.")
          message("Please type ?dtwclustControl for more information.\n")

          for (i in id_oldargs) {
               if (args[i] != "reps") {
                    var <- MYCALL[-1][[args[i]]]
                    slot(control, args[i]) <- if (is.numeric(var)) as.integer(var) else var

               } else
                    control@nrep <- as.integer(MYCALL$reps)

               ## remove from dots to avoid duplicate matches
               dots[[args[i]]] <- NULL
          }
     }

     ## ----------------------------------------------------------------------------------------------------------
     ## Preprocess
     ## ----------------------------------------------------------------------------------------------------------

     if (!is.null(preproc) && is.function(preproc)) {
          data <- preproc(data)
          preproc_char <- as.character(substitute(preproc))[[1]]

     } else if (type == "partitional" && is.character(centroid) && centroid == "shape") {
          preproc <- zscore
          preproc_char <- "zscore"
          data <- zscore(data)

     } else if (is.null(preproc)) {
          preproc <- function(x) x
          preproc_char <- "none"

     } else stop("Invalid preprocessing")

     ## ----------------------------------------------------------------------------------------------------------
     ## Further options
     ## ----------------------------------------------------------------------------------------------------------

     if (control@save.data)
          datalist <- data
     else
          datalist <- list()

     lengths <- sapply(data, length)
     diff_lengths <- length(unique(lengths)) > 1L

     consistency_check(distance, "dist", trace = control@trace, lengths = diff_lengths, silent = FALSE)

     if(type == "partitional" && diff_lengths)
          consistency_check(centroid, "cent", trace = control@trace)

     if (diff_lengths && distance %in% c("dtw", "dtw2"))
          control@symmetric <- FALSE
     else if (distance %in% c("lbk", "lbi"))
          control@symmetric <- FALSE
     else if (distance %in% c("sbd", "dtw", "dtw2"))
          control@symmetric <- TRUE

     ## For parallel computation
     control@packages <- c("dtwclust", control@packages)
     check_parallel() # register doSEQ if necessary

     if (type == "partitional") {

          ## =================================================================================================================
          ## Partitional
          ## =================================================================================================================

          if (k < 2L)
               stop("At least two clusters must be defined")

          if (is.character(centroid))
               centroid <- match.arg(centroid, c("mean", "median", "shape", "dba", "pam"))

          ## ----------------------------------------------------------------------------------------------------------
          ## Distance function
          ## ----------------------------------------------------------------------------------------------------------

          # dtwdistfun.R
          distfun <- dtwdistfun(distance = distance,
                                control = control,
                                distmat = distmat)

          ## ----------------------------------------------------------------------------------------------------------
          ## PAM precompute?
          ## ----------------------------------------------------------------------------------------------------------

          ## for a check near the end, changed if appropriate
          distmat_flag <- FALSE

          # precompute distance matrix?
          if (is.character(centroid) && centroid == "pam") {

               ## check if distmat was not provided and should be precomputed
               if (!is.null(distmat)) {
                    if ( nrow(distmat) != length(data) || ncol(distmat) != length(data) )
                         stop("Dimensions of provided cross-distance matrix don't correspond to length of provided data")

                    ## distmat was provided in call
                    distmat_flag <- TRUE

                    if (control@trace)
                         cat("\n\tDistance matrix provided...\n\n")

               } else if (control@pam.precompute) {
                    if (tolower(distance) == "dtw_lb")
                         warning("Using dtw_lb with control@pam.precompute = TRUE is not advised.")

                    if (control@trace)
                         cat("\n\tPrecomputing distance matrix...\n\n")

                    distmat <- do.call("distfun", c(list(data), dots))

                    ## Redefine new distmat
                    assign("distmat", distmat, envir = environment(distfun))

                    gc(FALSE)
               }

          } else {
               # replace any given distmat if centroid != "pam"
               distmat <- NULL
          }

          ## ----------------------------------------------------------------------------------------------------------
          ## Centroid function
          ## ----------------------------------------------------------------------------------------------------------

          ## Closures, all_cent.R
          allcent <- all_cent(case = centroid,
                              distmat = distmat,
                              distfun = distfun,
                              control = control)

          centroid <- as.character(substitute(centroid))[[1]]

          ## ----------------------------------------------------------------------------------------------------------
          ## Further options
          ## ----------------------------------------------------------------------------------------------------------

          family <- new("dtwclustFamily",
                        dist = distfun,
                        allcent = allcent,
                        preproc = preproc)

          if (control@trace && control@nrep > 1L)
               message("Tracing of repetitions might not be available if done in parallel.\n")

          if (!is.null(seed))
               set.seed(seed)

          ## ----------------------------------------------------------------------------------------------------------
          ## Cluster
          ## ----------------------------------------------------------------------------------------------------------

          if (control@nrep > 1L) {
               ## I need to re-register any custom distances in each parallel worker
               dist_entry <- proxy::pr_DB$get_entry(distance)

               export <- c("kcca.list", "consistency_check")

               kc.list <- foreach(i = 1:control@nrep,
                                  .combine = list,
                                  .multicombine = TRUE,
                                  .packages = control@packages,
                                  .export = export) %dorng% {
                                       if (!consistency_check(dist_entry$names[1], "dist"))
                                            do.call(proxy::pr_DB$set_entry, dist_entry)

                                       kc <- do.call("kcca.list",
                                                     c(dots,
                                                       list(x = data,
                                                            k = k,
                                                            family = family,
                                                            iter.max = control@iter.max,
                                                            trace = control@trace)))

                                       gc(FALSE)

                                       kc
                                  }
          } else {
               ## Just one repetition
               kc.list <- list(do.call("kcca.list",
                                       c(dots,
                                         list(x = data,
                                              k = k,
                                              family = family,
                                              iter.max = control@iter.max,
                                              trace = control@trace))))
          }

          ## ----------------------------------------------------------------------------------------------------------
          ## Prepare results
          ## ----------------------------------------------------------------------------------------------------------

          ## Replace distmat with NULL so that, if the distance function is called again, it won't subset it
          assign("distmat", NULL, envir = environment(family@dist))

          ## If distmat was provided, let it be shown in the results
          if (distmat_flag) {
               family@dist <- function(...) stop("'distmat' provided in call, no distance calculations performed")

               if (!is.null(attr(distmat, "method")))
                    distance <- attr(distmat, "method")
               else
                    distance <- "unknown"
          }

          ## Create objects
          dtwc <- lapply(kc.list, function(kc) {
               new("dtwclust",
                   call = MYCALL,
                   control = control,
                   family = family,
                   distmat = distmat,

                   type = type,
                   method = "NA",
                   distance = distance,
                   centroid = centroid,
                   preproc = preproc_char,

                   centers = kc$centers,
                   k = as.integer(k),
                   cluster = kc$cluster,

                   clusinfo = kc$clusinfo,
                   cldist = kc$cldist,
                   iter = kc$iter,
                   converged = kc$converged,

                   datalist = datalist)
          })

          if (control@nrep == 1L)
               RET <- dtwc[[1]]
          else
               RET <- dtwc

     } else if (type == "hierarchical") {

          ## =================================================================================================================
          ## Hierarchical
          ## =================================================================================================================

          hclust_methods <- match.arg(method,
                                      c("ward.D", "ward.D2", "single", "complete",
                                        "average", "mcquitty", "median", "centroid",
                                        "all"),
                                      several.ok = TRUE)

          if (any(hclust_methods == "all"))
               hclust_methods <- c("ward.D", "ward.D2", "single", "complete",
                                   "average", "mcquitty", "median", "centroid")

          if (!is.null(distmat)) {
               if ( nrow(distmat) != length(data) || ncol(distmat) != length(data) )
                    stop("Dimensions of provided cross-distance matrix don't correspond to length of provided data")

               if (control@trace)
                    cat("\n\tDistance matrix provided...\n")

               D <- distmat
               distfun <- function(...) stop("'distmat' provided in call, no distance calculations performed")

               if (!is.null(attr(distmat, "method")))
                    distance <- attr(distmat, "method")
               else
                    distance <- "unknown"

          } else {
               ## Take advantage of the function I defined for the partitional methods
               ## Which can do calculations in parallel if appropriate
               distfun <- dtwdistfun(distance = distance,
                                     control = control,
                                     distmat = NULL)

               if (control@trace)
                    cat("\n\tCalculating distance matrix...\n")

               ## single argument is to calculate whole distance matrix
               D <- do.call("distfun", c(list(data), dots))

          }

          ## Required form for 'hclust'
          Dist <- D[lower.tri(D)]

          ## Needed attribute for 'hclust' (case sensitive)
          attr(Dist, "Size") <- length(data)
          attr(Dist, "method") <- attr(D, "method")

          if (control@trace)
               cat("\n\tPerforming hierarchical clustering...\n\n")

          ## Cluster
          hc.list <- lapply(hclust_methods, function(method) {
               hc <- stats::hclust(Dist, method)

               ## cutree and corresponding centers
               cluster <- stats::cutree(hc, k)

               centers <- sapply(1:k, function(kcent) {
                    id_k <- cluster == kcent

                    d_sub <- D[id_k, id_k, drop = FALSE]

                    id_center <- which.min(apply(d_sub, 1, sum))

                    which(id_k)[id_center]
               })

               ## Some additional cluster information (taken from flexclust)
               cldist <- as.matrix(D[,centers][cbind(1:length(data), cluster)])
               size <- as.vector(table(cluster))
               clusinfo <- data.frame(size = size,
                                      av_dist = as.vector(tapply(cldist[,1], cluster, sum))/size)

               new("dtwclust", hc,
                   call = MYCALL,
                   control = control,
                   family = new("dtwclustFamily",
                                dist = distfun,
                                preproc = preproc),
                   distmat = D,

                   type = type,
                   method = method,
                   distance = distance,
                   centroid = "NA",
                   preproc = preproc_char,

                   centers = data[centers],
                   k = as.integer(k),
                   cluster = cluster,

                   clusinfo = clusinfo,
                   cldist = cldist,
                   iter = 1L,
                   converged = TRUE,

                   datalist = datalist)
          })

          if (length(hc.list) == 1L)
               RET <- hc.list[[1]]
          else
               RET <- hc.list

     } else if (type == "tadpole") {

          ## =================================================================================================================
          ## TADPole
          ## =================================================================================================================

          control@window.size <- consistency_check(control@window.size, "window")
          control@norm <- "L2"

          if (is.null(dc))
               stop("Please specify 'dc' for tadpole algorithm")
          if (dc < 0)
               stop("The cutoff distance 'dc' must be positive")

          ## ----------------------------------------------------------------------------------------------------------
          ## Check data
          ## ----------------------------------------------------------------------------------------------------------

          consistency_check(data, "tslist")

          ## ----------------------------------------------------------------------------------------------------------
          ## Cluster
          ## ----------------------------------------------------------------------------------------------------------

          if (control@trace)
               cat("\nEntering TADPole...\n")

          R <- TADPole(data, window.size = control@window.size, k = k, dc = dc, error.check = FALSE)

          if (control@trace) {
               cat("\nTADPole completed, pruning percentage = ",
                   formatC(100-R$distCalcPercentage, digits = 3, width = -1, format = "fg"),
                   "%\n\n",
                   sep = "")
          }

          ## ----------------------------------------------------------------------------------------------------------
          ## Prepare results
          ## ----------------------------------------------------------------------------------------------------------

          ## mainly for predict generic
          distfun <- dtwdistfun("dtw_lb", control = control, distmat = NULL)

          ## Some additional cluster information (taken from flexclust)
          subdistmat <- distfun(data, data[R$centers][R$cl], pairwise = TRUE)
          cldist <- as.matrix(subdistmat)
          size <- as.vector(table(R$cl))
          clusinfo <- data.frame(size = size,
                                 av_dist = as.vector(tapply(cldist[,1], R$cl, sum))/size)

          RET <- new("dtwclust",
                     call = MYCALL,
                     control = control,
                     family = new("dtwclustFamily",
                                  dist = distfun,
                                  preproc = preproc),
                     distmat = NULL,

                     type = type,
                     distance = "DTW2_LB",
                     centroid = "pam (TADPole)",
                     preproc = preproc_char,

                     centers = data[R$centers],
                     k = as.integer(k),
                     cluster = as.integer(R$cl),

                     clusinfo = clusinfo,
                     cldist = cldist,
                     iter = 1L,
                     converged = TRUE,

                     datalist = datalist)
     }

     toc <- proc.time() - tic

     if (class(RET) == "dtwclust")
          RET@proctime <- toc
     else
          RET <- lapply(RET, function(ret) {
               ret@proctime <- toc

               ret
          })

     if (type == "partitional" && control@nrep > 1L)
          attr(RET, "rng") <- attr(kc.list, "rng")

     if (control@trace)
          cat("\tElapsed time is", toc["elapsed"], "seconds.\n\n")

     RET
}
