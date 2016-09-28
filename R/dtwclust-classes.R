#' Class definition for \code{dtwclustFamily}
#'
#' Formal S4 class with a family of functions used in \code{\link{dtwclust}}.
#'
#' The custom implementations also handle parallelization.
#'
#' Since the distance function makes use of \code{proxy}, it also supports any extra \code{\link[proxy]{dist}}
#' parameters in \code{...}.
#'
#' The prototype includes the \code{cluster} function for partitional methods, as well as a
#' pass-through \code{preproc} function.
#'
#' @slot dist The function to calculate the distance matrices.
#' @slot allcent The function to calculate centroids on each iteration.
#' @slot cluster The function used to assign a series to a cluster.
#' @slot preproc The function used to preprocess the data (relevant for \code{\link[stats]{predict}}).
#'
#' @name dtwclustFamily-class
#' @rdname dtwclustFamily-class
#' @aliases dtwclustFamily
#'
#' @exportClass dtwclustFamily
#'
setClass("dtwclustFamily",

         slots = c(dist = "function",
                   allcent = "function",
                   cluster = "function",
                   preproc = "function"),

         prototype = prototype(preproc = function(x, ...) x,
                               cluster = function(distmat = NULL, ...) {
                                    if (is.null(distmat))
                                         stop("Something is wrong, couldn't calculate distances.")

                                    max.col(-distmat, "first")
                               })
)

#' Class definition for \code{dtwclustControl}
#'
#' Formal S4 class with several control parameters used in \code{\link{dtwclust}}.
#'
#' Default values are shown at the end.
#'
#' @section Common parameters:
#'
#' \itemize{
#'   \item \code{window.size} = \code{NULL}
#'   \item \code{norm} = "L1"
#'   \item \code{delta} = 1e-3
#'   \item \code{trace} = \code{FALSE}
#'   \item \code{save.data} = \code{TRUE}
#'   \item \code{symmetric} = \code{FALSE}
#'   \item \code{packages} = \code{character(0)}
#' }
#'
#' @slot window.size Integer or \code{NULL}. Window constraint for DTW, DBA and LB calculations.
#' @slot norm Character. Pointwise distance for DTW, DBA and the LBs. Either \code{"L1"} for Manhattan distance
#' or \code{"L2"} for Euclidean. Ignored for \code{distance = "DTW"} (which always uses \code{"L1"}) and
#' \code{distance = "DTW2"} (which always uses \code{"L2"}).
#' @slot delta Numeric. Stopping criterion for \code{\link{DBA}} centroids and fuzzy clustering.
#' @slot trace Logical flag. If \code{TRUE}, more output regarding the progress is printed to screen.
#' @slot save.data Return a "copy" of the data in the returned object? Because of the way \code{R} handles
#' things internally, all copies should point to the same memory address.
#' @slot symmetric Logical flag. Is the distance function symmetric? In other words, is \code{dist(x,y)} ==
#' \code{dist(y,x)}? If \code{TRUE}, only half the distance matrix needs to be computed.
#' @slot packages Character vector with the names of any packages required for custom \code{proxy} functions.
#' See Parallel Computing section in \code{\link{dtwclust}}.
#'
#' @section Only for partitional procedures:
#'
#' \itemize{
#'   \item \code{dba.iter} = \code{15L}
#'   \item \code{pam.precompute} = \code{TRUE}
#' }
#'
#' @slot dba.iter Integer. Maximum number of iterations for \code{\link{DBA}} centroids.
#' @slot pam.precompute Logical flag. Precompute the whole distance matrix once and reuse it on each iteration
#' if using PAM centroids. Otherwise calculate distances at every iteration.
#'
#' @section Only for fuzzy clustering:
#'
#' \itemize{
#'   \item \code{fuzziness} = \code{2}
#' }
#'
#' @slot fuzziness Numeric. Exponent used for fuzzy clustering. Commonly termed \code{m} in the literature.
#'
#' @section For both partitional and fuzzy:
#'
#' \itemize{
#'   \item \code{iter.max} = \code{100L}
#'   \item \code{nrep} = \code{1L}
#' }
#'
#' @slot iter.max Integer. Maximum number of iterations.
#' @slot nrep Integer. How many times to repeat clustering with different starting points. See section
#' Repetitions in \code{\link{dtwclust}}.
#'
#' @name dtwclustControl-class
#' @rdname dtwclustControl-class
#' @aliases dtwclustControl
#'
#' @exportClass dtwclustControl
#'
NULL

setClassUnion("intORnull", c("integer", "NULL"))

setClass("dtwclustControl",

         slots = c(window.size = "intORnull",
                   norm = "character",
                   dba.iter = "integer",
                   fuzziness = "numeric",
                   delta = "numeric",
                   pam.precompute = "logical",
                   iter.max = "integer",
                   trace = "logical",
                   nrep = "integer",
                   save.data = "logical",
                   symmetric = "logical",
                   packages = "character"),

         prototype = prototype(window.size = NULL,
                               norm = "L1",
                               dba.iter = 15L,
                               fuzziness = 2,
                               delta = 1e-3,
                               pam.precompute = TRUE,
                               iter.max = 100L,
                               trace = FALSE,
                               nrep = 1L,
                               save.data = TRUE,
                               symmetric = FALSE,
                               packages = character(0))
)

#' Class definition for \code{dtwclust}
#'
#' Formal S4 class.
#'
#' This class contains \code{\link[stats]{hclust}} as superclass and supports all its methods. Plot is a
#' special case (see \code{\link{dtwclust-methods}}).
#'
#' Please note that not all slots will contain valid information for all clustering types. In some cases,
#' for example for fuzzy and hierarchical clustering, some results are computed assuming a hard partition
#' is created based on the fuzzy memberships or dendrogram tree, and the provided value of \code{k}.
#'
#' @slot call The function call.
#' @slot control An object of class \code{\link{dtwclustControl}}.
#' @slot family An object of class \code{\link{dtwclustFamily}}.
#' @slot distmat If computed, the cross-distance matrix.
#' @slot k Integer indicating the number of desired clusters.
#' @slot cluster Integer vector indicating which cluster a series belongs to (crisp partition).
#' @slot fcluster Numeric matrix that contains membership of fuzzy clusters. It has one row for each
#' series and one column for each cluster. The rows must sum to 1. Only relevant for fuzzy clustering.
#' @slot iter The number of iterations used.
#' @slot converged A logical indicating whether the function converged.
#' @slot clusinfo A data frame with two columns: \code{size} indicates the number of series each cluster has,
#' and \code{av_dist} indicates, for each cluster, the average distance between series and their respective
#' centroids (crisp partition).
#' @slot centers Deprecated, please use \code{centroids} instead.
#' @slot centroids A list with the centroid time series.
#' @slot cldist A column vector with the distance between each series in the data and its corresponding
#' centroid (crisp partition).
#' @slot type A string indicating one of the supported clustering types of \code{\link{dtwclust}}.
#' @slot method A string indicating which hierarchical method was used.
#' @slot distance A string indicating the distance used.
#' @slot centroid A string indicating the centroid used.
#' @slot preproc A string indicating the preprocessing used.
#' @slot datalist The provided data in the form of a list, where each element is a time series.
#' @slot proctime Time during function execution, as measured with \code{\link[base]{proc.time}}.
#' @slot dots The contents of the original call's ellipsis (...).
#'
#' @name dtwclust-class
#' @rdname dtwclust-class
#'
#' @exportClass dtwclust
#'
NULL

setClass("proc_time4", contains = "numeric", slots = c(names = "character"))
setOldClass("proc_time", S4Class = "proc_time4")
removeClass("proc_time4")

setClass("hclust4", contains = "list", slots = c(names = "character"))
setOldClass("hclust", S4Class = "hclust4")
removeClass("hclust4")

setClass("dtwclust", contains = c("hclust"),
         slots = c(call = "call",
                   control = "dtwclustControl",
                   family = "dtwclustFamily",
                   distmat = "ANY",

                   k = "integer",
                   cluster = "integer",
                   fcluster = "matrix",
                   iter = "integer",
                   converged = "logical",
                   clusinfo = "data.frame",

                   centers = "list",
                   centroids = "list",
                   cldist = "matrix",

                   type = "character",
                   method = "character",
                   distance = "character",
                   centroid = "character",
                   preproc = "character",
                   datalist = "list",
                   proctime = "proc_time",
                   dots = "list"))
