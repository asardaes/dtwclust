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
#' TODO
#'
setClass("tsclustFamily",
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

#' Class definition for \code{TSClusters} and derived classes
#'
#' Formal S4 classes for time-series clusters.
#'
#' @rdname TSClusters-class
#' @exportClass TSClusters
#'
#' @details
#'
#' TBD
#'
setClass("TSClusters",
         slots = c(call = "call",
                   family = "tsclustFamily",
                   control = "list",
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
                   dots = "list"))

#' @rdname TSClusters-class
#' @exportClass PartitionalTSClusters
#'
#' @details
#'
#' TBD2
#'
setClass("PartitionalTSClusters", contains = c("TSClusters"),
         slots = c(iter = "integer",
                   converged = "logical",
                   clusinfo = "data.frame",
                   cldist = "matrix"))

#' @rdname TSClusters-class
#' @exportClass HierarchicalTSClusters
#'
#' @details
#'
#' TBD3
#'
setClass("HierarchicalTSClusters", contains = c("TSClusters", "hclust"),
         slots = c(method = "character",
                   clusinfo = "data.frame",
                   cldist = "matrix"))

#' @rdname TSClusters-class
#' @exportClass FuzzyTSClusters
#'
#' @details
#'
#' TBD4
#'
setClass("FuzzyTSClusters", contains = c("TSClusters"),
         slots = c(iter = "integer",
                   converged = "logical",
                   fcluster = "matrix"))
