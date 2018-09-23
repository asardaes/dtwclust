#' @section Proxy version:
#'
#'   The version registered with \code{\link[proxy]{dist}} is custom (\code{loop = FALSE} in
#'   \code{\link[proxy]{pr_DB}}). The custom function handles multi-threaded parallelization
#'   directly (with \code{\link[RcppParallel:RcppParallel-package]{RcppParallel}}). It uses all
#'   available threads by default (see
#'   \code{\link[RcppParallel:defaultNumThreads]{RcppParallel::defaultNumThreads()}}), but this can
#'   be changed by the user with
#'   \code{\link[RcppParallel:setThreadOptions]{RcppParallel::setThreadOptions()}}.
#'
#'   An exception to the above is when it is called within a \code{\link[foreach:foreach]{foreach}}
#'   parallel loop \strong{made by dtwclust}. If the parallel workers do not have the number of
#'   threads explicitly specified, this function will default to 1 thread per worker. See the
#'   parallelization vignette for more information (\code{browseVignettes("dtwclust")}).
#'
