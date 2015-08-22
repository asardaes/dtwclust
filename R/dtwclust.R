#' Time series clustering under DTW
#'
#' This function uses the DTW distance and related lower bounds to cluster time series. For now, all series
#' must be univariate and, in the case of partitional methods, have equal lengths.
#'
#' Partitional algorithms are implemented via \code{\link[flexclust]{kcca}}. Hierarchical algorithms use the
#' \code{\link[stats]{hclust}} function. The \code{tadpole} algorithm uses the \code{\link{TADPole}} function.
#'
#' In case of partitional algorithms, \code{data} should be in the form of a matrix. In the other cases,
#' it may be a matrix or a list, but the matrix will be coerced to a list. A matrix input requires that all
#' time series have equal lengths. If the lengths vary slightly between time series, reinterpolating them to
#' a common length is most likely an acceptable approach (Ratanamahatana and Keogh, 2004). If this is not the
#' case, then clustering them directly is probably ill-advised. See the examples.
#'
#' @section Distance:
#'
#' If a custom distance function is provided, it will receive the data as the first argument. For partitional
#' algorithms, the second argument will be the cluster centers (i.e. other time series) in the form of a matrix
#' where each row is a center series. If hierarchical algorithms are used, the function will also receive the
#' elements of \code{...}.
#'
#' For partitional algorithms, the function \emph{could} make use of the \code{window.size} and \code{norm}
#' parameters, which \emph{should} be detected thanks to \code{R}'s lexical scoping, however this cannot
#' be guaranteed.
#'
#' The function should return a distance matrix, ideally of class \code{crossdist}. In case of partitional
#' algorithms, the time series in the data should be along the rows, and the cluster centers along the
#' columns of the distance matrix.
#'
#' The other option is to provide a string. For partitional algorithms, it can be a supported distance of
#' \code{\link[flexclust]{kccaFamily}}. For hierarchical, it can be any distance function registered in
#' \code{\link[proxy]{dist}}. In the latter case, all extra parameters should be provided in \code{...}.
#'
#' Additionally, with either type of algorithm, it can be one of the following
#' custom implementations:
#'
#' \itemize{
#'   \item \code{"dtw"}: DTW with L1 norm and optionally a Sakoe-Chiba constraint.
#'   \item \code{"dtw2"}: DTW with L2 norm and optionally a Sakoe-Chiba constraint.
#'   \item \code{"dtw_lb"}: DTW with L1 or L2 norm and optionally a Sakoe-Chiba constraint. Some computations
#'   are avoided by first estimating the distance matrix with Lemire's lower bound and then iteratively
#'   refining with DTW. See \code{\link{dtw_lb}}.
#'   \item \code{"lbk"}: Keogh's lower bound with either L1 or L2 norm for the Sakoe-Chiba constraint.
#'   \item \code{"lbi"}: Lemire's lower bound with either L1 or L2 norm for the Sakoe-Chiba constraint.
#'   \item \code{"sbd"}: Shape-based distance. Each series is z-normalized in this case. As a result,
#'   the cluster centers (for partitional methods) are also z-normalized. See \code{\link{SBD}} for more
#'   details.
#' }
#'
#' @section Centroid:
#'
#' In the case of partitional algorithms, a suitable function should calculate the cluster centers. In this
#' case, the centers are themselves time series.
#'
#' If a custom function is provided, it will receive a matrix as only argument. Each row will be a time series
#' that belongs to a given cluster. The function should return a numeric vector with the center time series.
#'
#' The other option is to provide a character string. The following options are available:
#'
#' \itemize{
#'   \item \code{"mean"}: The average along each dimension. In other words, the average of all \eqn{x^j_i}
#'   among the \eqn{j} series that belong to the same cluster for all time points \eqn{t_i}.
#'   \item \code{"median"}: The median along each dimension. Similar to mean.
#'   \item \code{"shape"}: Shape averaging. See \code{\link{shape_extraction}} for more details.
#'   \item \code{"dba"}: DTW Barycenter Averaging. See \code{\link{DBA}} for more details.
#'   \item \code{"pam"}: Partition around medoids. This basically means that the cluster centers are always
#'   one of the time series in the data. In this case, the distance matrix is pre-computed once using all
#'   time series in the data and then re-used at each iteration.
#' }
#'
#' @section Sakoe-Chiba Constraint:
#'
#' A global constraint to speed up the DTW calculation is the Sakoe-Chiba band (Sakoe and Chiba, 1978). To
#' use it, a window width must be defined via \code{window.size}.
#'
#' The windowing constraint uses a centered window. The calculations expect a value in \code{window.size}
#' that represents the distance between the point considered and one of the edges of the window. Therefore,
#' if, for example, \code{window.size = 10}, the warping for an observation \eqn{x_i} considers the points
#' between \eqn{x_{i-10}} and \eqn{x_{i+10}}, resulting in \code{10*2 + 1 = 21} observations falling within
#' the window.
#'
#' @section Preprocessing:
#'
#' It is strongly advised to use z-normalization in case of \code{centroid = "shape"}, because the resulting
#' series have this normalization (see \code{\link{shape_extraction}}). The user can, however, specify a
#' custom function that performs any transformation on the data, but the user must make sure that the format
#' stays consistent, i.e. that a matrix where each row is a series for partitional methods, or a list of time
#' series for hierarchical methods. For example, the z-normalization could be implemented as
#' \code{t(apply(data, 1, zscore))} or \code{lapply(data, zscore)} respectively.
#'
#' The function will receive the data as first argument and, in case hierarchical methods are used, the
#' contents of \code{...} as the second argument.
#'
#' @section Notes:
#'
#' In order to ensure that the parameter values are detected correctly by the included functions
#' when partitional clustering is used, the
#' environment of the \code{dtwclust} function is assigned as an attribute of \code{data} via
#' \code{attr(data, "env")} \code{<-} \code{environment()}. If the user alters the dataset with a
#' preprocessing function, it should make sure that this attribute is maintained.
#'
#' Notice that the lower bounds are defined only for time series of equal lengths. \code{DTW} and \code{DTW2}
#' don't require this, but they are much slower to compute.
#'
#' The lower bounds are \strong{not} symmetrical, and \code{DTW} is only symmetrical if series are of equal
#' lengths.
#'
#' Specifying \code{distance = "sbd"} and \code{centroid = "shape"} is equivalent to the k-Shape algorithm
#' (Papparizos and Gravano, 2015). See \code{\link{SBD}} and \code{\link{shape_extraction}} for more info.
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
#' @examples
#'
#' # Load data
#' data(uciCT)
#'
#' # Reinterpolate to same length and coerce as matrix
#' data <- t(sapply(CharTraj, reinterpolate, newLength = 205))
#'
#' # Simple partitional clustering with L2 distance and PAM
#' kc.l2 <- dtwclust(data, k = 20, distance = "kmeans", centroid = "pam",
#'                   seed = 3247, trace = TRUE)
#' cat("Rand index for L2+PAM:", randIndex(kc.l2, CharTrajLabels), "\n\n")
#'
#' # TADPole clustering (takes around 5 seconds)
#' kc.tadp <- dtwclust(data, type = "tadpole", k = 20,
#'                     window.size = 20, dc = 1.5, save.data = TRUE)
#' cat("Rand index for TADPole:", randIndex(kc.tadp, CharTrajLabels), "\n\n")
#' plot(kc.tadp)
#'
#' # Modify plot
#' plot(kc.tadp, cl = 1:4, labs.arg = list(title = "TADPole, clusters 1 through 4",
#'                                         x = "time", y = "series"))
#'
#' \dontrun{
#' # Hierarchical clustering based on shabe-based distance
#' hc.sbd <- dtwclust(data, type = "hierarchical", distance = "sbd")
#' cl.sbd <- cutree(hc.sbd, 20)
#' cat("Rand index for HC+SBD:", randIndex(cl.sbd, CharTrajLabels), "\n\n")
#'
#' # Use full DTW and PAM (takes around two minutes)
#' kc.dtw <- dtwclust(data, k = 20, seed = 3251, trace = TRUE)
#'
#' # Use full DTW with DBA centroids (takes around five minutes)
#' kc.dba <- dtwclust(data, k = 20, centroid = "dba", seed = 3251, trace = TRUE)
#' }
#'
#' @author Alexis Sarda-Espinosa
#'
#' @param data A list where each element is a time series, or a numerical matrix where each row is a time
#' series. Only the latter is supported in case of \code{type = "partitional"}.
#' @param type What type of clustering method to use, \code{partitional}, \code{hierarchical} or \code{tadpole}.
#' @param k Numer of desired clusters in partitional methods.
#' @param method Which linkage method to use in hierarchical methods. See \code{\link[stats]{hclust}}.
#' @param distance One of the supported distance definitions (see Distance section). Ignored for
#' \code{type = "tadpole"}.
#' @param centroid Either a supported string or an appropriate function to calculate centroids
#' when using partitional methods (see Centroid section).
#' @param preproc Function to preprocess data. Defaults to \code{zscore} \emph{only} if \code{centroid =
#' "shape"}, but will be replaced by a custom function if provided. See Preprocessing section.
#' @param window.size Window constraint for DTW and LB calculations. See Sakoe-Chiba section.
#' @param norm Pointwise distance for DTW and LB. Either \code{L1} for Manhattan distance or \code{L2}
#' for Euclidean. Ignored for \code{distance = "DTW"} (which always uses \code{L1}) and
#' \code{distance = "DTW2"} (which always uses \code{L2}).
#' @param dc Cutoff distance for TADPole algorithm.
#' @param dba.iter Maximum number of iterations for \code{\link{DBA}} centroids.
#' @param control Parameters for partitional clustering algorithms. See
#' \code{\link[flexclust]{flexclustControl}}.
#' @param save.data Return a copy of the data in the returned object? Ignored for hierarchical clustering.
#' @param seed Random seed for reproducibility of partitional algorithms.
#' @param trace Boolean flag. If true, more output regarding the progress is printed to screen.
#' @param ... Additional arguments to pass to \code{\link[proxy]{dist}} or a custom function.
#'
#' @return An object with formal class \code{\link{dtwclust-class}} if \code{type = "partitional" | "tadpole"}. Otherwise
#' an object with class \code{hclust} as returned by \code{\link[stats]{hclust}}.
#'
#' @export
#' @import flexclust
#' @importFrom stats median
#' @importFrom stats hclust
#' @importFrom proxy dist
#' @importFrom modeltools ModelEnvMatrix

dtwclust <- function(data = NULL, type = "partitional", k = 2, method = "average",
                     distance = "dtw", centroid = "pam", preproc = NULL,
                     window.size = NULL, norm = "L1", dc = NULL,
                     dba.iter = 50, control = NULL, save.data = FALSE,
                     seed = NULL, trace = FALSE,
                     ...)
{
     ## =================================================================================================================
     ## Start
     ## =================================================================================================================

     if (is.null(data))
          stop("No data provided")

     tic <- proc.time()

     type <- match.arg(type, c("partitional", "hierarchical", "tadpole"))
     norm <- match.arg(norm, c("L1", "L2"))

     if (type == "partitional") {

          ## =================================================================================================================
          ## Partitional
          ## =================================================================================================================

          if (k < 2)
               stop("At least two clusters must be defined")

          ## Used by some of the custom functions so that they know where to look for parameters
          ## This is done automatically due to lexical scoping, but I rather do it explicitly
          attr(data, "env") <- environment()

          if (is.function(centroid)) {
               cent <- centroid
          } else {
               centroid <- match.arg(centroid, c("mean", "median", "shape", "dba", "pam"))

               cent <- switch(EXPR = centroid,

                              mean = function(x) { apply(x, 2, mean) },

                              median = function(x) { apply(x, 2, stats::median) },

                              shape = shape_extraction, # placeholder, will be ignored (see utils.R)

                              dba = DBA, # placeholder, will be ignored (see utils.R)

                              pam = function(x) NULL, # placeholder, will be ignored (see utils.R)

                              ## Otherwise
                              stop("Unsupported centroid calculation method"))
          }

          if (is.function(distance)) {
               family <- kccaFamily(name = as.character(substitute(distance)),
                                    dist = distance,
                                    cent = cent)

          } else if (is.character(distance)) {
               family <- kccaFamilies(distance, cent, window.size, norm)

          } else {
               stop("Unspported distance definition")
          }

          ## ----------------------------------------------------------------------------------------------------------
          ## Replace specific functions if necessary
          ## ----------------------------------------------------------------------------------------------------------

          distmat <- NULL

          if (centroid == "shape") {
               family@allcent <- allcent_se

               family@preproc <- preproc_se

          } else if (centroid == "dba") {
               family@allcent <- allcent_dba

          } else if (centroid == "pam") {
               family@allcent <- allcent_pam

               distmat <- preproc_pam(data, family)

          }

          if (!is.null(preproc)) {
               if (is.function(preproc))
                    family@preproc <- preproc
               else
                    stop("Invalid preprocessing")
          }

          ## ----------------------------------------------------------------------------------------------------------
          ## Further options
          ## ----------------------------------------------------------------------------------------------------------

          if (is.null(control))
               ctrl <- new("flexclustControl")
          else
               ctrl <- as(control, "flexclustControl")

          if (trace)
               ctrl@verbose <- 1

          if (!is.null(seed))
               set.seed(seed)

          ## ----------------------------------------------------------------------------------------------------------
          ## Cluster
          ## ----------------------------------------------------------------------------------------------------------

          kc <- flexclust::kcca(x = data,
                                k = k,
                                family = family,
                                simple = TRUE,
                                control = ctrl,
                                save.data = save.data)

          ## ----------------------------------------------------------------------------------------------------------
          ## Prepare results
          ## ----------------------------------------------------------------------------------------------------------

          toc <- proc.time() - tic

          if (trace)
               cat("\n\tElapsed time is", toc["elapsed"], "seconds.\n\n")

          dtwc <- new("dtwclust", kc,
                      type = type,
                      distance = ifelse(is.function(distance), as.character(substitute(distance)), distance),
                      centroid = ifelse(is.function(centroid), as.character(substitute(centroid)), centroid))

          dtwc

     } else if (type == "hierarchical") {

          ## =================================================================================================================
          ## Hierarchical
          ## =================================================================================================================

          if (trace)
               cat("\n\tCalculating distance matrix...\n")

          if (is.matrix(data))
               x <- lapply(seq_len(nrow(data)), function(i) data[i,])
          else if (is.list(data))
               x <- data
          else
               stop("Unsupported type for data")

          if (!is.null(preproc) && is.function(preproc)) {
               data <- preproc(data, ...)
          }

          if (is.function(distance)) {
               D <- distance(data, ...)

          } else if (is.character(distance)) {
               D <- switch(EXPR = distance,

                           ## Full DTW with L1 norm
                           dtw = {
                                if (is.null(window.size)) {
                                     d <- proxy::dist(x = x, y = x,
                                                      method = "DTW", dist.method = "L1",
                                                      ...)
                                } else {
                                     d <- proxy::dist(x = x, y = x,
                                                      method = "DTW", dist.method = "L1",
                                                      window.type = "slantedband", window.size = window.size,
                                                      ...)
                                }

                                d
                           },

                           ## Full DTW with L2 norm
                           dtw2 = {
                                if (is.null(window.size)) {
                                     d <- proxy::dist(x = x, y = x,
                                                      method = "DTW2",
                                                      ...)
                                } else {
                                     d <- proxy::dist(x = x, y = x,
                                                      method = "DTW2",
                                                      window.type = "slantedband", window.size = window.size,
                                                      ...)
                                }

                                d
                           },

                           ## DTW with aid of lower bounds
                           dtw_lb = {
                                window.size <- consistency_check(window.size, "window")

                                dtw_lb(x, x, window.size, norm = norm, error.check=TRUE)
                           },


                           ## Lemire's improved lower bound with L1
                           lbi = {
                                window.size <- consistency_check(window.size, "window")

                                proxy::dist(x = x, y = x,
                                            method = "LBI", window.size = window.size, norm = norm,
                                            force.symmetry = TRUE, error.check=TRUE,
                                            ...)
                           },

                           ## Keogh's lower bound but with L1
                           lbk = {
                                window.size <- consistency_check(window.size, "window")

                                proxy::dist(x = x, y = x,
                                            method = "LBK", window.size = window.size, norm = norm,
                                            force.symmetry = TRUE, error.check=TRUE,
                                            ...)
                           },

                           ## Paparrizos' shape-based distance
                           sbd = {
                                x <- lapply(x, zscore)

                                proxy::dist(x = x, y = x,
                                            method = "SBD",
                                            ...)
                           },


                           ## Otherwise
                           proxy::dist(x = x, y = x, method = distance, ...)
               )

          } else {
               stop("Unspported distance definition")
          }

          if (trace)
               cat("\n\tPerforming hierarchical clustering...\n")

          ## Required form for 'hclust'
          D <- D[lower.tri(D)]

          ## Needed attribute for 'hclust' (case sensitive)
          attr(D, "Size") <- length(x)

          hc <- stats::hclust(D, method = method)

          toc <- proc.time() - tic

          if (trace)
               cat("\n\tElapsed time is", toc["elapsed"], "seconds.\n\n")

          hc

     } else if (type == "tadpole") {

          ## =================================================================================================================
          ## TADPole
          ## =================================================================================================================

          window.size <- consistency_check(window.size, "window")

          if (is.null(dc))
               stop("The user must specify 'dc' for this method")
          if (dc < 0)
               stop("The cutoff distance 'dc' must be positive")

          ## ----------------------------------------------------------------------------------------------------------
          ## Adjust inputs
          ## ----------------------------------------------------------------------------------------------------------

          if (is.matrix(data))
               x <- lapply(seq_len(nrow(data)), function(i) data[i,])
          else if (is.list(data))
               x <- data
          else
               stop("Unsupported type for data")

          consistency_check(x, "tslist")

          ## ----------------------------------------------------------------------------------------------------------
          ## Cluster
          ## ----------------------------------------------------------------------------------------------------------

          if (trace)
               cat("\nEntering TADPole...\n")

          R <- TADPole(x, window.size = window.size, k = k, dc = dc, error.check = FALSE)

          ## ----------------------------------------------------------------------------------------------------------
          ## Prepare results
          ## ----------------------------------------------------------------------------------------------------------

          if (save.data) {
               tadpc <- new("dtwclust",
                            type = type,
                            distance = "dtw",
                            centroid = "TADPole (PAM)",

                            centers = data[R$centers, ],
                            k = as.integer(k),
                            cluster = as.integer(R$cl),
                            data = modeltools::ModelEnvMatrix(designMatrix = data))
          } else {
               tadpc <- new("dtwclust",
                            type = type,
                            distance = "dtw",
                            centroid = "TADPole (PAM)",

                            centers = data[R$centers, ],
                            k = as.integer(k),
                            cluster = as.integer(R$cl))
          }

          toc <- proc.time() - tic

          if (trace)
               cat("\n\tElapsed time is", toc["elapsed"], "seconds.\n\n")

          tadpc
     }
}
