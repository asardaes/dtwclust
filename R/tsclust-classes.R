#' Class definition for \code{tsclustFamily}
#'
#' Formal S4 class with a family of functions used in \code{\link{tsclust}}.
#'
#' @exportClass tsclustFamily
#'
#' @details
#'
#' The custom implementations also handle parallelization.
#'
#' Since the distance function makes use of \code{proxy}, it also supports any extra
#' \code{\link[proxy]{dist}} parameters in \code{...}.
#'
#' The prototype includes the \code{cluster} function for partitional methods, as well as a
#' pass-through \code{preproc} function.
#'
#' @slot dist The function to calculate the distance matrices.
#' @slot allcent The function to calculate centroids on each iteration.
#' @slot cluster The function used to assign a series to a cluster.
#' @slot preproc The function used to preprocess the data (relevant for
#'   \code{\link[stats]{predict}}).
#'
#' @examples
#'
#' # The dist() function in tsclustFamily works like proxy::dist() but supports
#' # parallelization and optimized symmetric calculations. If you like, you can
#' # use the function more or less directly, but provide a control argument when
#' # creating the family.
#'
#' \dontrun{
#' data(uciCT)
#' fam <- new("tsclustFamily", dist = "gak",
#'            control = partitional_control(symmetric = TRUE))
#'
#' crossdist <- fam@dist(CharTraj, window.size = 18L)
#' }
#'
#' # If you want the fuzzy family, use fuzzy = TRUE
#' ffam <- new("tsclustFamily", control = fuzzy_control(), fuzzy = TRUE)
#'
#'
setClass("tsclustFamily",
         slots = c(dist = "function",
                   allcent = "function",
                   cluster = "function",
                   preproc = "function"),
         prototype = prototype(preproc = function(x, ...) { x },

                               cluster = function(distmat = NULL, ...) {
                                   max.col(-distmat, "first")
                               })
)

#' Class definition for \code{TSClusters} and derived classes
#'
#' Formal S4 classes for time-series clusters. See class hierarchy and slot organization at the
#' bottom.
#'
#' @rdname TSClusters-class
#' @exportClass TSClusters
#'
#' @details
#'
#' The base class is \code{TSClusters}. The 3 classes that inherit from it are:
#' \code{PartitionalTSClusters}, \code{HierarchicalTSClusters} and \code{FuzzyTSClusters}.
#'
#' \code{HierarchicalTSClusters} also contain \code{\link[stats]{hclust}} as parent class.
#'
#' Package \pkg{clue} is supported, but generics from \pkg{flexclust} are not. See also
#' \code{\link{tsclusters-methods}}.
#'
#' If you want to transform a \code{\link{dtwclust-class}} object to TSClusters, just use
#' \code{as(dtwclust_obj, "TSClusters")}, although it may not work perfectly.
#'
#' @slot call The function call.
#' @slot family An object of class \code{\link{tsclustFamily-class}}.
#' @slot control An appropriate control object for \code{\link{tsclust}}. See
#'   \code{\link{tsclust-controls}}.
#' @slot datalist The provided data in the form of a list, where each element is a time series.
#' @slot type A string indicating one of the supported clustering types of \code{\link{tsclust}}.
#' @slot distance A string indicating the distance used.
#' @slot centroid A string indicating the centroid used.
#' @slot preproc A string indicating the preprocessing used.
#' @slot k Integer indicating the number of desired clusters.
#' @slot cluster Integer vector indicating which cluster a series belongs to (crisp partition). For
#'   fuzzy clustering, this is based on \strong{distance}, not on \code{fcluster}. For hierarchical,
#'   this is obtained by calling \code{\link[stats]{cutree}} with the given value of \code{k}.
#' @slot centroids A list with the centroid time series.
#' @slot distmat If computed, the cross-distance matrix.
#' @slot proctime Time during function execution, as measured with \code{\link[base]{proc.time}}.
#' @slot dots The contents of the original call's ellipsis (...).
#' @slot args The contents of the original call's \code{args} parameter. See
#'   \code{\link{tsclust_args}}.
#'
#' @section TSClusters:
#'
#'   The base class contains the following slots:
#'
#'   \itemize{
#'     \item \code{call}
#'     \item \code{family}
#'     \item \code{control}
#'     \item \code{datalist}
#'     \item \code{type}
#'     \item \code{distance}
#'     \item \code{centroid}
#'     \item \code{preproc}
#'     \item \code{k}
#'     \item \code{cluster}
#'     \item \code{centroids}
#'     \item \code{distmat}
#'     \item \code{proctime}
#'     \item \code{dots}
#'     \item \code{args}
#'   }
#'
#' @seealso
#'
#' \code{\link{tsclusters-methods}}
#'
setClass("TSClusters",
         slots = c(call = "call",
                   family = "tsclustFamily",
                   control = "ANY",
                   datalist = "list",

                   type = "character",
                   distance = "character",
                   centroid = "character",
                   preproc = "character",

                   k = "integer",
                   cluster = "integer",
                   centroids = "list",
                   distmat = "ANY",

                   proctime = "proc_time",
                   dots = "list",
                   args = "ANY"))

#' @rdname TSClusters-class
#' @exportClass PartitionalTSClusters
#'
#' @slot iter The number of iterations used.
#' @slot converged A logical indicating whether the function converged.
#' @slot clusinfo A data frame with two columns: \code{size} indicates the number of series each
#'   cluster has, and \code{av_dist} indicates, for each cluster, the average distance between
#'   series and their respective centroids (crisp partition).
#' @slot cldist A column vector with the distance between each series in the data and its
#'   corresponding centroid (crisp partition).
#'
#' @section PartitionalTSClusters:
#'
#'   This class adds the following slots to the base class:
#'
#'   \itemize{
#'     \item \code{iter}
#'     \item \code{converged}
#'     \item \code{clusinfo}
#'     \item \code{cldist}
#'   }
#'
setClass("PartitionalTSClusters", contains = c("TSClusters"),
         slots = c(iter = "integer",
                   converged = "logical",
                   clusinfo = "data.frame",
                   cldist = "matrix"))

#' @rdname TSClusters-class
#' @exportClass HierarchicalTSClusters
#'
#' @slot method A string indicating which hierarchical method was used.
#'
#' @section HierarchicalTSClusters:
#'
#'   This class adds the following slots to the base class:
#'
#'   \itemize{
#'     \item \code{method}
#'     \item \code{clusinfo}
#'     \item \code{cldist}
#'   }
#'
setClass("HierarchicalTSClusters", contains = c("TSClusters", "hclust"),
         slots = c(method = "character",
                   clusinfo = "data.frame",
                   cldist = "matrix"))

#' @rdname TSClusters-class
#' @exportClass FuzzyTSClusters
#'
#' @slot fcluster Numeric matrix that contains membership of fuzzy clusters. It has one row for each
#'   series and one column for each cluster. The rows must sum to 1. Only relevant for fuzzy
#'   clustering.
#'
#' @section FuzzyTSClusters:
#'
#'   This class adds the following slots to the base class:
#'
#'   \itemize{
#'     \item \code{iter}
#'     \item \code{converged}
#'     \item \code{fcluster}
#'   }
#'
setClass("FuzzyTSClusters", contains = c("TSClusters"),
         slots = c(iter = "integer",
                   converged = "logical",
                   fcluster = "matrix"))
