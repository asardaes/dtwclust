#' Time series clustering
#'
#' This is the main function to perform time series clustering. It supports partitional, hierarchical, fuzzy
#' and TADPole clustering. See the details and the examples for more information.
#'
#' Partitional and fuzzy clustering procedures use a custom implementation. Hierarchical clustering is done
#' with \code{\link[stats]{hclust}}. TADPole clustering uses the \code{\link{TADPole}} function. Specifying
#' \code{type} = \code{"partitional"}, \code{distance} = \code{"sbd"} and \code{centroid} = \code{"shape"} is
#' equivalent to the k-Shape algorithm (Paparrizos and Gravano, 2015).
#'
#' The \code{data} may be a matrix, a data frame or a list. Matrices and data frames are coerced to a list,
#' the former row-wise and the latter column-wise. Only lists can have series with different lengths. Most of
#' the optimizations require series to have the same length, so consider reinterpolating them to save some
#' time (see Ratanamahatana and Keogh, 2004; \code{\link{reinterpolate}}). No missing values are allowed.
#'
#' Several parameters can be adjusted with the \code{control} argument. See \code{\link{dtwclustControl}}. In
#' the following sections, elements marked with an asterisk (*) are those that can be adjutsed with this
#' argument.
#'
#' @section Partitional Clustering:
#'
#' Stochastic algorithm that creates a hard partition of the data into \code{k} clusters, where each cluster
#' has a center/centroid. In case of time series clustering, the centroids are also time series.
#'
#' The cluster centroids are first randomly initialized by selecting some of the series in the data. The
#' distance between each series and each centroid is calculated, and the series are assigned to the cluster
#' whose centroid is closest. The centroids are then updated by using a given rule, and the procedure is
#' repeated until no series changes from one cluster to another, or until the maximum number of iterations*
#' has been reached. The distance and centroid definitions can be specified through the corresponding
#' parameters of this function. See their respective sections below.
#'
#' Note that it is possible for a cluster to become empty, in which case a new cluster is reinitialized
#' randomly. However, if you see that the algorithm doesn't converge or the overall distance sum increases,
#' it could mean that the chosen value of \code{k} is too large, or the chosen distance measure is not
#' able to assess similarity effectively. The random reinitialization attempts to enforce a certain number
#' of clusters, but can result in instability in the aforementioned cases.
#'
#' @section Fuzzy Clustering:
#'
#' This procedure is very similar to partitional clustering, except that each series no longer belongs
#' exclusively to one cluster, but belongs to each cluster to a certain degree. For each series, the total
#' degree of membership across clusters must sum to 1.
#'
#' The default implementation uses the fuzzy c-means algorithm. In its definition, an objective function is to
#' be minimized. The objective is defined in terms of a squared distance, which is usually the Euclidean (L2)
#' distance, although the definition could be modified. The \code{distance} parameter of this function controls
#' the one that is utilized. The fuzziness of the clustering can be controlled by means of the fuzziness
#' exponent*. Bear in mind that the centroid definition of fuzzy c-means requires equal dimensions, which
#' means that all series must have the same length. This problem can be circumvented by applying
#' transformations to the series (see for example D'Urso and Maharaj, 2009).
#'
#' Note that the fuzzy clustering could be transformed to a crisp one by finding the highest membership
#' coefficient. Some of the slots of the object returned by this function assume this, so be careful with
#' interpretation (see \code{\link{dtwclust-class}}).
#'
#' @section Hierarchical Clustering:
#'
#' This is a deterministic algorithm that creates a hierarchy of groups by using different linkage methods
#' (see \code{\link[stats]{hclust}}). The linkage method is controlled through the \code{method} parameter
#' of this function, with the additional option "all" that uses all of the available methods. The distance
#' to be used can be controlled with the \code{distance} parameter.
#'
#' The hierarchy does not imply a specific number of clusters, but one can be induced by cutting the resulting
#' dendrogram (see \code{\link[stats]{cutree}}). As with fuzzy clustering, this results in a crisp partition,
#' and some of the slots of the returned object are calculated by cutting the dendrogram so that \code{k}
#' clusters are created.
#'
#' @section TADPole Clustering:
#'
#' TADPole clustering adopts a relatively new clustering framework and adapts it to time series clustering
#' with DTW. Because of the way it works, it can be considered a kind of Partitioning Around Medoids (PAM).
#' This means that the cluster centers are always elements of the data. However, this algorithm is
#' deterministic, depending on the value of the cutoff distance \code{dc}, which can be controlled with the
#' corresponding parameter of this function.
#'
#' The algorithm relies on the DTW lower bounds, which are only defined for time series of equal length.
#' Additionally, it requires a window constraint* for DTW. See the Sakoe-Chiba constraint section below.
#' Unlike the other algorithms, TADPole always uses DTW2 as a distance (with a symmetric1 step pattern).
#'
#' @section Centroid Calculation:
#'
#' In the case of partitional/fuzzy algorithms, a suitable function should calculate the cluster centers at
#' every iteration. In this case, the centers are themselves time series. Fuzzy clustering uses the standard
#' fuzzy c-means centroid by default.
#'
#' In either case, a custom function can be provided. If one is provided, it will receive the following
#' inputs in the shown order (examples for partitional clustering are shown in parenthesis):
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
#' Therefore, the function should \emph{always} include the ellipsis \code{...} in its definition. In case of
#' fuzzy clustering, the membership vectors (2nd and 5th elements above) are matrices with number of rows
#' equal to amount of elements in the data, and number of columns equal to the number of desired clusters.
#' Each row must sum to 1.
#'
#' The other option is to provide a character string for the custom implementations. The following options
#' are available:
#'
#' \itemize{
#'   \item "mean": The average along each dimension. In other words, the average of all \eqn{x^j_i}
#'   among the \eqn{j} series that belong to the same cluster for all time points \eqn{t_i}.
#'   \item "median": The median along each dimension. Similar to mean.
#'   \item "shape": Shape averaging. By default, all series are z-normalized in this case, since the resulting
#'   centroids will also have this normalization. See \code{\link{shape_extraction}} for more details.
#'   \item "dba": DTW Barycenter Averaging. See \code{\link{DBA}} for more details.
#'   \item "pam": Partition around medoids (PAM). This basically means that the cluster centers are always
#'   one of the time series in the data. In this case, the distance matrix can be pre-computed once using all
#'   time series in the data and then re-used at each iteration. It usually saves overhead overall.
#'   \item "fcm": Fuzzy c-means. Only supported for fuzzy clustering and always used for that type of clustering
#'   if a string is provided in \code{centroid}.
#' }
#'
#' These check for the special cases where parallelization might be desired. Note that only \code{shape},
#' \code{dba} and \code{pam} support series of different length. Also note that, for \code{shape} and
#' \code{dba}, this support has a caveat: the final centroids' length will depend on the length of those
#' series that were randomly chosen at the beginning of the clustering algorithm. For example, if the series
#' in the dataset have a length of either 10 or 15, 2 clusters are desired, and the initial choice selects
#' two series with length of 10, the final centroids will have this same length.
#'
#' @section Distance Measures:
#'
#' The distance measure to be used with partitional, hierarchical and fuzzy clustering can be modified with
#' the \code{distance} parameter. The supported option is to provide a string, which must represent a
#' compatible distance registered with \code{proxy}'s \code{\link[proxy]{dist}}. Registration is done via
#' \code{\link[proxy]{pr_DB}}, and extra parameters can be provided in \code{...}.
#'
#' Note that you are free to create your own distance functions and register them. Optionally, you can use
#' one of the following custom implementations (all registered with \code{proxy}):
#'
#' \itemize{
#'   \item \code{"dtw"}: DTW with L1 norm and optionally a Sakoe-Chiba/Slanted-band constraint*.
#'   \item \code{"dtw2"}: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint*.
#'   \item \code{"dtw_lb"}: DTW with L1 or L2 norm* and optionally a Sakoe-Chiba constraint*. Some
#'   computations are avoided by first estimating the distance matrix with Lemire's lower bound and then
#'   iteratively refining with DTW. See \code{\link{dtw_lb}}. Not suitable for \code{pam.precompute}* =
#'   \code{TRUE}.
#'   \item \code{"lbk"}: Keogh's lower bound with either L1 or L2 norm* for the Sakoe-Chiba constraint*.
#'   \item \code{"lbi"}: Lemire's lower bound with either L1 or L2 norm* for the Sakoe-Chiba constraint*.
#'   \item \code{"sbd"}: Shape-based distance. See \code{\link{SBD}} for more details.
#' }
#'
#' DTW2 is done with \code{\link[dtw]{dtw}}, but it differs from the result you would obtain if you specify
#' \code{L2} as \code{dist.method}: with \code{DTW2}, pointwise distances (the local cost matrix) are
#' calculated with \code{L1} norm, \emph{each} element of the matrix is squared and the result is fed into
#' \code{\link[dtw]{dtw}}, which finds the optimum warping path. The square root of the resulting
#' distance is \emph{then} computed.
#'
#' Only \code{dtw}, \code{dtw2} and \code{sbd} support series of different length. The lower bounds
#' are probably unsuitable for direct clustering unless series are very easily distinguishable.
#'
#' If you create your own distance, register it with \code{proxy}, and it includes the ellipsis (\code{...})
#' in its definition, it will receive the following parameters*:
#'
#' \itemize{
#'   \item \code{window.type}: Either \code{"none"} for a \code{NULL} \code{window.size}, or \code{"slantedband"}
#'   otherwise
#'   \item \code{window.size}: The provided window size
#'   \item \code{norm}: The provided desired norm
#'   \item \code{...}: Any additional parameters provided in the original call's ellipsis
#' }
#'
#' Whether your function makes use of them or not, is up to you.
#'
#' If you know that the distance function is symmetric, and you use a hierarchical algorithm or a partitional
#' algorithm with PAM centroids and \code{pam.precompute}* = \code{TRUE}, some time can be saved by
#' calculating only half the distance matrix. Therefore, consider setting the symmetric* control parameter
#' to \code{TRUE} if this is the case.
#'
#' @section Sakoe-Chiba Constraint:
#'
#' A global constraint to speed up the DTW calculations is the Sakoe-Chiba band (Sakoe and Chiba, 1978). To
#' use it, a window width* must be defined.
#'
#' The windowing constraint uses a centered window. The function expects a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10(2) + 1 = 21} observations falling within
#' the window.
#'
#' The computations actually use a \code{slantedband} window, which is equivalent to the Sakoe-Chiba one
#' if series have equal length, and stays along the diagonal of the local cost matrix if series have
#' different length.
#'
#' @section Preprocessing:
#'
#' It is strongly advised to use z-normalization in case of \code{centroid = "shape"}, because the resulting
#' series have this normalization (see \code{\link{shape_extraction}}). Therefore, \code{\link{zscore}} is the
#' default in this case. The user can, however, specify a custom function that performs any transformation
#' on the data, but the user must make sure that the format stays consistent, i.e. a list of time series.
#'
#' Setting to \code{NULL} means no preprocessing (except for \code{centroid = "shape"}). A provided function
#' will receive the data as first and only argument.
#'
#' It is convenient to provide this function if you're planning on using the \code{\link[stats]{predict}}
#' generic.
#'
#' @section Repetitions:
#'
#' Due to their stochastic nature, partitional clustering is usually repeated* several times with different
#' random seeds to allow for different starting points. This function uses \code{\link[rngtools]{RNGseq}} to
#' obtain different seed streams for each repetition, utilizing the \code{seed} parameter (if provided) to
#' initialize it. If more than one repetition is made, the streams are returned in an attribute called
#' \code{rng}.
#'
#' Technically, you can also perform random repetitions for fuzzy clustering, although it might be difficult
#' to evaluate the results, since they are usually evaluated relative to each other and not in an absolute
#' way. Ideally, the groups wouldn't change too much once the algorithm converges.
#'
#' Multiple values of \code{k} can also be provided to get different partitions using any \code{type} of
#' clustering.
#'
#' Repetitions are greatly optimized when PAM centroids are used and the whole distance matrix is precomputed*,
#' since said matrix is reused for every repetition, and can be comptued in parallel (see Parallel section).
#'
#' @section Parallel Computing:
#'
#' Please note that running tasks in parallel does \strong{not} guarantee faster computations.
#' The overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#' The user can register a parallel backend, for eample with the \code{doParallel} package, in order to do
#' the repetitions in parallel, as well as distance and some centroid calculations.
#'
#' Unless each repetition requires a few seconds, parallel computing probably isn't worth it. As such, I
#' would only use this feature with \code{shape} and \code{DBA} centroids, or an expensive distance function
#' like \code{DTW}.
#'
#' If you register a parallel backend, the function will also try to do the calculation of the distance
#' matrices in parallel. This should work with any function registered with \code{\link[proxy]{dist}} via
#' \code{\link[proxy]{pr_DB}} whose \code{loop} flag is set to \code{TRUE}. If the function requires special
#' packages to be loaded, provide their names in the \code{packages}* slot of \code{control}. Note that
#' "dtwclust" is always loaded in each parallel worker, so that doesn't need to be included. Alternatively,
#' you may want to pre-load \code{dtwclust} in each worker with \code{\link[parallel]{clusterEvalQ}}.
#'
#' In case of multiple repetitions, each worker gets a repetition task. Otherwise, the tasks (which can be
#' a distance matrix or a centroid calculation) are usually divided into chunks according to the number of
#' workers available.
#'
#' @section Notes:
#'
#' The lower bounds are defined only for time series of equal length. \code{DTW} and \code{DTW2}
#' don't require this, but they are much slower to compute.
#'
#' The lower bounds are \strong{not} symmetric, and \code{DTW} is only symmetric if series have equal
#' lengths.
#'
#' @references
#'
#' Sakoe H and Chiba S (1978). ``Dynamic programming algorithm optimization for spoken word
#' recognition.'' \emph{Acoustics, Speech and Signal Processing, IEEE Transactions on},
#' \strong{26}(1), pp. 43-49. ISSN 0096-3518, \url{http://doi.org/10.1109/TASSP.1978.1163055}.
#'
#' Ratanamahatana A and Keogh E (2004). ``Everything you know about dynamic time warping is wrong.''
#' In \emph{3rd Workshop on Mining Temporal and Sequential Data, in conjunction with 10th ACM SIGKDD
#' Int. Conf. Knowledge Discovery and Data Mining (KDD-2004), Seattle, WA}.
#'
#' Bedzek, J.C. (1981). Pattern recognition with fuzzy objective function algorithms.
#'
#' D'Urso, P., & Maharaj, E. A. (2009). Autocorrelation-based fuzzy clustering of time series.
#' Fuzzy Sets and Systems, 160(24), 3565-3589.
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In \emph{Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data},
#' series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' @example inst/dtwclust-examples.R
#'
#' @seealso
#'
#' \code{\link{dtwclust-methods}}, \code{\link{dtwclust-class}}, \code{\link{dtwclustControl}},
#' \code{\link{dtwclustFamily}}.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @param data A list where each element is a time series, or a numeric matrix (it will be coerced to a list
#' row-wise).
#' @param type What type of clustering method to use: \code{partitional}, \code{hierarchical}, \code{tadpole}
#' or \code{fuzzy}.
#' @param k Numer of desired clusters. It may be a numeric vector with different values.
#' @param method One or more linkage methods to use in hierarchical procedures. See \code{\link[stats]{hclust}}.
#' You can provide a character vector to compute different hierarchical cluster structures in one go, or
#' specify \code{method} \code{=} \code{"all"} to use all the available ones.
#' @param distance A supported distance from \code{proxy}'s \code{\link[proxy]{dist}} (see Distance section).
#' Ignored for \code{type} = \code{"tadpole"}.
#' @param centroid Either a supported string or an appropriate function to calculate centroids
#' when using partitional methods. See Centroid section.
#' @param preproc Function to preprocess data. Defaults to \code{zscore} \emph{only} if \code{centroid}
#' \code{=} \code{"shape"}, but will be replaced by a custom function if provided. See Preprocessing section.
#' @param dc Cutoff distance for the \code{\link{TADPole}} algorithm.
#' @param control Named list of parameters or \code{dtwclustControl} object for clustering algorithms. See
#' \code{\link{dtwclustControl}}. \code{NULL} means defaults.
#' @param seed Random seed for reproducibility of partitional and fuzzy algorithms.
#' @param distmat If a cross-distance matrix is already available, it can be provided here so it's re-used.
#' Only relevant if \code{centroid} = "pam" or \code{type} = "hierarchical". See examples.
#' @param ... Additional arguments to pass to \code{\link[proxy]{dist}} or a custom function.
#'
#' @return An object with formal class \code{\link{dtwclust-class}}.
#'
#' If \code{control@nrep > 1} and a partitional procedure is used, \code{length(method)} \code{> 1} and
#' hierarchical procedures are used, or \code{length(k)} \code{>} \code{1}, a list of objects is returned.
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

     type <- match.arg(type, c("partitional", "hierarchical", "tadpole", "fuzzy"))

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

     args <- names(MYCALL[-1L])
     id_oldargs <- which(args %in% slotNames("dtwclustControl"))
     if (length(id_oldargs) > 0L) {
          message("Control arguments should not be provided in ... to prevent duplicate matches.")
          message("Please type ?dtwclustControl for more information.\n")

          for (i in id_oldargs) {
               var <- MYCALL[-1L][[args[i]]]
               slot(control, args[i]) <- var

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

     Lengths <- lengths(data)
     diff_lengths <- length(unique(Lengths)) > 1L

     consistency_check(distance, "dist", trace = control@trace, Lengths = diff_lengths, silent = FALSE)

     if(type %in% c("partitional", "fuzzy")) {
          if (is.character(centroid)) {
               if (type == "fuzzy")
                    centroid <- "fcm"
               else
                    centroid <- match.arg(centroid, c("mean", "median", "shape", "dba", "pam"))
          }

          if (diff_lengths)
               consistency_check(centroid, "cent", trace = control@trace)
     }

     if (diff_lengths && distance %in% c("dtw", "dtw2"))
          control@symmetric <- FALSE
     else if (distance %in% c("lbk", "lbi"))
          control@symmetric <- FALSE
     else if (distance %in% c("sbd", "dtw", "dtw2"))
          control@symmetric <- TRUE

     ## For parallel computation
     control@packages <- c("dtwclust", control@packages)
     check_parallel() # register doSEQ if necessary

     if (any(k < 2L))
          stop("At least two clusters must be defined")
     if (any(k > length(data)))
          stop("Cannot have more clusters than series in the data")

     if (type %in% c("partitional", "fuzzy")) {

          ## =================================================================================================================
          ## Partitional or fuzzy
          ## =================================================================================================================

          ## ----------------------------------------------------------------------------------------------------------
          ## Distance function
          ## ----------------------------------------------------------------------------------------------------------

          # ddist.R
          distfun <- ddist(distance = distance,
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
                              control = control,
                              fuzzy = isTRUE(type == "fuzzy"))

          centroid <- as.character(substitute(centroid))[[1L]]

          ## ----------------------------------------------------------------------------------------------------------
          ## Further options
          ## ----------------------------------------------------------------------------------------------------------

          family <- new("dtwclustFamily",
                        dist = distfun,
                        allcent = allcent,
                        preproc = preproc)

          if (type == "fuzzy")
               family@cluster <- fcm_cluster # utils.R

          if (control@trace && control@nrep > 1L)
               message("Tracing of repetitions might not be available if done in parallel.\n")

          ## ----------------------------------------------------------------------------------------------------------
          ## Cluster
          ## ----------------------------------------------------------------------------------------------------------

          if (length(k) == 1L && control@nrep == 1L) {
               rngtools::setRNG(rngtools::RNGseq(1L, seed = seed, simplify = TRUE))

               ## Just one repetition
               kc.list <- list(do.call("kcca.list",
                                       c(dots,
                                         list(x = data,
                                              k = k,
                                              family = family,
                                              control = control,
                                              fuzzy = isTRUE(type == "fuzzy")))))

          } else {
               ## I need to re-register any custom distances in each parallel worker
               dist_entry <- proxy::pr_DB$get_entry(distance)

               export <- c("kcca.list", "consistency_check")

               rng <- rngtools::RNGseq(length(k) * control@nrep, seed = seed, simplify = FALSE)
               rng <- lapply(parallel::splitIndices(length(rng), length(k)), function(i) rng[i])

               comb0 <- if(control@nrep > 1L) c else list

               k0 <- k
               rng0 <- rng

               kc.list <- foreach(k = k0, rng = rng0, .combine = comb0, .multicombine = TRUE,
                                  .packages = control@packages, .export = export) %:%
                    foreach(i = 1L:control@nrep, .combine = list, .multicombine = TRUE,
                            .packages = control@packages, .export = export) %dopar% {
                                 rngtools::setRNG(rng[[i]])

                                 if (!consistency_check(dist_entry$names[1], "dist"))
                                      do.call(proxy::pr_DB$set_entry, dist_entry)

                                 kc <- do.call("kcca.list",
                                               c(dots,
                                                 list(x = data,
                                                      k = k,
                                                      family = family,
                                                      control = control,
                                                      fuzzy = isTRUE(type == "fuzzy"))))

                                 gc(FALSE)

                                 kc
                            }
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
          RET <- lapply(kc.list, function(kc) {
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
                   k = kc$k,
                   cluster = kc$cluster,
                   fcluster = kc$fcluster,

                   clusinfo = kc$clusinfo,
                   cldist = kc$cldist,
                   iter = kc$iter,
                   converged = kc$converged,

                   datalist = datalist)
          })

          if (length(RET) == 1L) RET <- RET[[1L]]

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
               distfun <- ddist(distance = distance,
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
          hc <- lapply(hclust_methods, function(method) stats::hclust(Dist, method))

          RET <- lapply(k, function(k) {
               lapply(hc, function(hc) {
                    ## cutree and corresponding centers
                    cluster <- stats::cutree(hc, k)

                    centers <- sapply(1L:k, function(kcent) {
                         id_k <- cluster == kcent

                         d_sub <- D[id_k, id_k, drop = FALSE]

                         id_center <- which.min(apply(d_sub, 1L, sum))

                         which(id_k)[id_center]
                    })

                    ## Some additional cluster information (taken from flexclust)
                    cldist <- as.matrix(D[,centers][cbind(1L:length(data), cluster)])
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
                        method = hc$method,
                        distance = distance,
                        centroid = "NA",
                        preproc = preproc_char,

                        centers = data[centers],
                        k = as.integer(k),
                        cluster = cluster,
                        fcluster = matrix(NA_real_),

                        clusinfo = clusinfo,
                        cldist = cldist,
                        iter = 1L,
                        converged = TRUE,

                        datalist = datalist)
               })
          })

          RET <- unlist(RET, recursive = FALSE)

          if (length(RET) == 1L) RET <- RET[[1L]]

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

          ## mainly for predict generic
          distfun <- ddist("dtw_lb", control = control, distmat = NULL)

          RET <- foreach(k = k, .combine = list, .multicombine = TRUE, .packages = "dtwclust") %dopar% {
               R <- TADPole(data, window.size = control@window.size, k = k, dc = dc, error.check = FALSE)

               if (control@trace) {
                    cat("\nTADPole completed, pruning percentage = ",
                        formatC(100 - R$distCalcPercentage, digits = 3L, width = -1L, format = "fg"),
                        "%\n\n",
                        sep = "")
               }

               ## ----------------------------------------------------------------------------------------------------------
               ## Prepare results
               ## ----------------------------------------------------------------------------------------------------------

               ## Some additional cluster information (taken from flexclust)
               subdistmat <- distfun(data, data[R$centers][R$cl], pairwise = TRUE)
               cldist <- as.matrix(subdistmat)
               size <- as.vector(table(R$cl))
               clusinfo <- data.frame(size = size,
                                      av_dist = as.vector(tapply(cldist[ , 1L], R$cl, sum))/size)

               new("dtwclust",
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
                   fcluster = matrix(NA_real_),

                   clusinfo = clusinfo,
                   cldist = cldist,
                   iter = 1L,
                   converged = TRUE,

                   datalist = datalist)
          }
     }

     toc <- proc.time() - tic

     if (class(RET) == "dtwclust")
          RET@proctime <- toc
     else
          RET <- lapply(RET, function(ret) {
               ret@proctime <- toc

               ret
          })

     if (type %in% c("partitional", "fuzzy") && (control@nrep > 1L || length(k) > 1L))
          attr(RET, "rng") <- unlist(rng0, recursive = FALSE, use.names = FALSE)

     if (control@trace)
          cat("\tElapsed time is", toc["elapsed"], "seconds.\n\n")

     RET
}
