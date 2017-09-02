#' Class definition for `tsclustFamily`
#'
#' Formal S4 class with a family of functions used in [tsclust()].
#'
#' @exportClass tsclustFamily
#'
#' @details
#'
#' The custom implementations also handle parallelization.
#'
#' Since the distance function makes use of `proxy`, it also supports any extra [proxy::dist()]
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
#'   The dist() function in here works like [proxy::dist()] but supports parallelization and
#'   optimized symmetric calculations. If you like, you can use the function more or less directly,
#'   but provide a control argument when creating the family. However, bear in mind the following
#'   considerations.
#'
#'   - The second argument is called `centroids` (inconsistent with [proxy::dist()]).
#'   - If `control$distmat` is *not* `NULL`, the function will try to subset it.
#'   - If `control$symmetric` is `TRUE`, `centroids` is `NULL`, *and* there is no argument
#'     `pairwise` that is `TRUE`, only half the distance matrix will be computed.
#'     + If the distance was registered in [proxy::pr_DB] with `loop = TRUE` and more than one
#'       parallel worker is detected, the computation will be in parallel, otherwise it will be
#'       sequential with [proxy::dist()].
#'   - The function always returns a `crossdist` matrix.
#'
#'   Note that all distances implemented as part of \pkg{dtwclust} have custom proxy loops, so see
#'   their respective documentation to see what optimizations apply to each one.
#'
#'   For distances *not* included in \pkg{dtwclust}, the symmetric, parallel case mentioned above
#'   makes chunks for parallel workers, but they are not perfectly balanced, so some workers might
#'   finish before the others.
#'
#' @section Centroid function:
#'
#'   The default partitional allcent() function is a closure with the implementations of the
#'   included centroids. The ones for [DBA()] and [shape_extraction()] can use parallelization. Its
#'   formal arguments are described in the Centroid Calculation section from [tsclust()].
#'
#' @note
#'
#' This class is meant to group together the relevant functions, but they are **not** linked with
#' each other automatically. In other words, neither `dist` nor `allcent` apply `preproc`. They
#' essentially don't know of each other's existence.
#'
#' @examples
#'
#' \dontrun{
#' data(uciCT)
#' # See "GAK" documentation
#' fam <- new("tsclustFamily", dist = "gak")
#'
#' # This is done with symmetric optimizations, regardless of control$symmetric
#' crossdist <- fam@dist(CharTraj, window.size = 18L)
#'
#' # This is done without symmetric optimizations, regardless of control$symmetric
#' crossdist <- fam@dist(CharTraj, CharTraj, window.size = 18L)
#'
#' # For non-dtwclust distances, symmetric optimizations only apply
#' # with an appropriate control AND a single data argument:
#' fam <- new("tsclustFamily", dist = "dtw",
#'            control = partitional_control(symmetric = TRUE))
#' fam@dist(CharTraj[1L:5L])
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

# ==================================================================================================
# Custom initialize
# ==================================================================================================

setMethod("initialize", "tsclustFamily",
          function(.Object, dist, allcent, ..., control = list(), fuzzy = FALSE) {
              dots <- list(...)
              dots$.Object <- .Object

              if (!missing(dist)) {
                  if (is.character(dist))
                      dots$dist <- ddist2(dist, control)
                  else
                      dots$dist <- dist
              }

              if (fuzzy) {
                  dots$cluster <- fcm_cluster # fuzzy.R
                  if (!missing(allcent) && is.character(allcent))
                      allcent <- match.arg(allcent, c("fcm", "fcmdd"))
              }

              if (!missing(allcent)) {
                  if (is.character(allcent))
                      dots$allcent <- all_cent2(allcent, control)
                  else if (is.function(allcent))
                      dots$allcent <- allcent
                  else
                      stop("Centroid definition must be either a function or a character")
              }

              do.call(methods::callNextMethod, dots, TRUE)
          })
