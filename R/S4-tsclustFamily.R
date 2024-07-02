#' Class definition for `tsclustFamily`
#'
#' Formal S4 class with a family of functions used in [tsclust()].
#'
#' @exportClass tsclustFamily
#' @importFrom methods setClass
#'
#' @details
#'
#' The custom implementations also handle parallelization.
#'
#' Since the distance function makes use of \pkg{proxy}, it also supports any extra [proxy::dist()]
#' parameters in `...`.
#'
#' The prototype includes the `cluster` function for partitional methods, as well as a pass-through
#' `preproc` function. The initializer expects a control from [tsclust-controls]. See more below.
#'
#' @slot dist The function to calculate the distance matrices.
#' @slot allcent The function to calculate centroids on each iteration.
#' @slot cluster The function used to assign a series to a cluster.
#' @slot preproc The function used to preprocess the data (relevant for [stats::predict()]).
#'
#' @section Distance function:
#'
#'   The family's dist() function works like [proxy::dist()] but supports parallelization and
#'   optimized symmetric calculations. If you like, you can use the function more or less directly,
#'   but provide a control argument when creating the family (see examples). However, bear in mind
#'   the following considerations.
#'
#'   - The second argument is called `centroids` (inconsistent with [proxy::dist()]).
#'   - If `control$distmat` is *not* `NULL`, the function will try to subset it.
#'   - If `control$symmetric` is `TRUE`, `centroids` is `NULL`, *and* there is no argument
#'     `pairwise` that is `TRUE`, only half the distance matrix will be computed.
#'
#'   Note that all distances implemented as part of \pkg{dtwclust} have custom proxy loops that use
#'   multi-threading independently of \pkg{foreach}, so see their respective documentation to see
#'   what optimizations apply to each one.
#'
#'   For distances *not* included in \pkg{dtwclust}, the computation can be in parallel using
#'   multi-processing with [foreach::foreach()]. If you install and load or attach (see
#'   [base::library()] or [base::loadNamespace()]) the \pkg{bigmemory} package, the function will
#'   take advantage of said package when all of the following conditions are met, reducing the
#'   overhead of data copying across processes:
#'
#'   - `control$symmetric` is `TRUE`
#'   - `centroids` is `NULL`
#'   - `pairwise` is `FALSE` or `NULL`
#'   - The distance was registered in [proxy::pr_DB] with `loop = TRUE`
#'   - A parallel backend with more than 1 worker has been registered with \pkg{foreach}
#'
#'   This symmetric, parallel case makes chunks for parallel workers, but they are not perfectly
#'   balanced, so some workers might finish before the others.
#'
#' @section Centroid function:
#'
#'   The default partitional allcent() function is a closure with the implementations of the
#'   included centroids. The ones for [DBA()], [shape_extraction()] and [sdtw_cent()] can use
#'   multi-process parallelization with [foreach::foreach()]. Its formal arguments are described in
#'   the Centroid Calculation section from [tsclust()].
#'
#' @note
#'
#' This class is meant to group together the relevant functions, but they are **not** linked with
#' each other automatically. In other words, neither `dist` nor `allcent` apply `preproc`. They
#' essentially don't know of each other's existence.
#'
#' @seealso
#'
#' [dtw_basic()], [dtw_lb()], [gak()], [lb_improved()], [lb_keogh()], [sbd()], [sdtw()].
#'
#' @examples
#'
#' \dontrun{
#' data(uciCT)
#' # See "GAK" documentation
#' fam <- new("tsclustFamily", dist = "gak")
#'
#' # This is done with symmetric optimizations, regardless of control$symmetric
#' crossdist <- fam@@dist(CharTraj, window.size = 18L)
#'
#' # This is done without symmetric optimizations, regardless of control$symmetric
#' crossdist <- fam@@dist(CharTraj, CharTraj, window.size = 18L)
#'
#' # For non-dtwclust distances, symmetric optimizations only apply
#' # with an appropriate control AND a single data argument:
#' fam <- new("tsclustFamily", dist = "dtw",
#'            control = partitional_control(symmetric = TRUE))
#' fam@@dist(CharTraj[1L:5L])
#'
#' # If you want the fuzzy family, use fuzzy = TRUE
#' ffam <- new("tsclustFamily", control = fuzzy_control(), fuzzy = TRUE)
#' }
#'
tsclustFamily <- methods::setClass(
    "tsclustFamily",
    slots = c(dist = "function",
              allcent = "function",
              cluster = "function",
              preproc = "function"),
    prototype = prototype(preproc = function(x, ...) { x },
                          cluster = function(distmat = NULL, ...) {
                              max.col(-distmat, "first")
                          })
)

# ==================================================================================================
# Membership update for fuzzy clustering
# ==================================================================================================

f_cluster <- function(distmat, m) {
    cprime <- apply(distmat, 1L, function(dist_row) { sum( (1 / dist_row) ^ (2 / (m - 1)) ) })
    u <- 1 / apply(distmat, 2L, function(dist_col) { cprime * dist_col ^ (2 / (m - 1)) })
    if (is.null(dim(u))) u <- rbind(u) # for predict generic
    u[is.nan(u)] <- 1 # in case fcmdd is used
    u
}

# ==================================================================================================
# Custom initialize
# ==================================================================================================

#' @importFrom methods as
#' @importFrom methods callNextMethod
#' @importFrom methods initialize
#' @importFrom methods setMethod
#' @importFrom rlang enexprs
#' @importFrom rlang env_bind
#'
setMethod("initialize", "tsclustFamily",
          function(.Object, dist, allcent, ..., control = list(), fuzzy = FALSE) {
              rlang::env_bind(environment(), ...)
              dots <- rlang::enexprs(...)
              dots$.Object <- quote(.Object)
              if (!missing(dist)) {
                  if (is.character(dist)) dist <- ddist2(dist, control)
                  dots$dist <- quote(dist)
              }
              if (fuzzy) {
                  dots$cluster <- quote(f_cluster)
                  if (!missing(allcent) && is.character(allcent))
                      allcent <- match.arg(allcent, c("fcm", "fcmdd"))
              }
              if (!missing(allcent)) {
                  if (is.character(allcent)) {
                      if (allcent %in% c("pam", "fcmdd")) {
                          if (!is.null(control$distmat) && !inherits(control$distmat, "Distmat")) {
                              control$distmat <- methods::as(control$distmat, "Distmat")
                          }
                      }
                      allcent <- all_cent2(allcent, control)
                  }
                  else if (!is.function(allcent)) {
                      stop("Centroid definition must be either a function or a character")
                  }
                  dots$allcent <- quote(allcent)
              }
              do.call(methods::callNextMethod, dots)
          })
