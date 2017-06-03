#' Class definition for `TSClusters` and derived classes
#'
#' Formal S4 classes for time-series clusters. See class hierarchy and slot organization at the
#' bottom.
#'
#' @rdname TSClusters-class
#' @exportClass TSClusters
#' @include tsclust-family.R
#'
#' @details
#'
#' The base class is `TSClusters`. The 3 classes that inherit from it are: `PartitionalTSClusters`,
#' `HierarchicalTSClusters` and `FuzzyTSClusters`.
#'
#' `HierarchicalTSClusters` also contain [stats::hclust()] as parent class.
#'
#' Package \pkg{clue} is supported, but generics from \pkg{flexclust} are not. See also
#' [tsclusters-methods].
#'
#' If you want to transform a [dtwclust-class] object to TSClusters, just use:
#'
#' `as(dtwclust_obj, "TSClusters")`
#'
#' although it may not work perfectly.
#'
#' @slot call The function call.
#' @slot family An object of class [tsclustFamily-class].
#' @slot control An appropriate control object for [tsclust()]. See [tsclust-controls].
#' @slot datalist The provided data in the form of a list, where each element is a time series.
#' @slot type A string indicating one of the supported clustering types of [tsclust()].
#' @slot distance A string indicating the distance used.
#' @slot centroid A string indicating the centroid used.
#' @slot preproc A string indicating the preprocessing used.
#' @slot k Integer indicating the number of desired clusters.
#' @slot cluster Integer vector indicating which cluster a series belongs to (crisp partition). For
#'   fuzzy clustering, this is based on **distance**, not on `fcluster`. For hierarchical, this is
#'   obtained by calling [stats::cutree()] with the given value of `k`.
#' @slot centroids A list with the centroid time series.
#' @slot distmat If computed, the cross-distance matrix.
#' @slot proctime Time during function execution, as measured with [base::proc.time()].
#' @slot dots The contents of the original call's ellipsis (...).
#' @slot args The contents of the original call's `args` parameter. See [tsclust_args()].
#' @slot seed The random seed that was used.
#'
#' @section TSClusters:
#'
#'   The base class contains the following slots:
#'
#'   - `call`
#'   - `family`
#'   - `control`
#'   - `datalist`
#'   - `type`
#'   - `distance`
#'   - `centroid`
#'   - `preproc`
#'   - `k`
#'   - `cluster`
#'   - `centroids`
#'   - `distmat`
#'   - `proctime`
#'   - `dots`
#'   - `args`
#'   - `seed`
#'
#' @seealso
#'
#' [tsclusters-methods]
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
                   args = "ANY",
                   seed = "integer"))

#' @rdname TSClusters-class
#' @exportClass PartitionalTSClusters
#'
#' @slot iter The number of iterations used.
#' @slot converged A logical indicating whether the function converged.
#' @slot clusinfo A data frame with two columns: `size` indicates the number of series each cluster
#'   has, and `av_dist` indicates, for each cluster, the average distance between series and their
#'   respective centroids (crisp partition).
#' @slot cldist A column vector with the distance between each series in the data and its
#'   corresponding centroid (crisp partition).
#'
#' @section PartitionalTSClusters:
#'
#'   This class adds the following slots to the base class:
#'
#'   - `iter`
#'   - `converged`
#'   - `clusinfo`
#'   - `cldist`
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
#'   - `method`
#'   - `clusinfo`
#'   - `cldist`
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
#'   - `iter`
#'   - `converged`
#'   - `fcluster`
#'
setClass("FuzzyTSClusters", contains = c("TSClusters"),
         slots = c(iter = "integer",
                   converged = "logical",
                   fcluster = "matrix"))
