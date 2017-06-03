## Needed old classes
methods::setClass("proc_time4", contains = "numeric", slots = c(names = "character"))
methods::setOldClass("proc_time", S4Class = "proc_time4")
methods::removeClass("proc_time4")

methods::setClass("hclust4", contains = "list", slots = c(names = "character"))
methods::setOldClass("hclust", S4Class = "hclust4")
methods::removeClass("hclust4")

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
#' `preproc` function.
#'
#' @slot dist The function to calculate the distance matrices.
#' @slot allcent The function to calculate centroids on each iteration.
#' @slot cluster The function used to assign a series to a cluster.
#' @slot preproc The function used to preprocess the data (relevant for [stats::predict()]).
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

              do.call(methods::callNextMethod, dots)
          })
