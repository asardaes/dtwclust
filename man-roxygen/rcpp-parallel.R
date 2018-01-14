#' @section Parallel Computing:
#'
#'   Please note that running tasks in parallel does \strong{not} guarantee faster computations. The
#'   overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#'   This function uses the \code{\link[RcppParallel:RcppParallel-package]{RcppParallel}} package
#'   for parallelization. It uses all available threads by default (see
#'   \code{\link[RcppParallel:defaultNumThreads]{RcppParallel::defaultNumThreads()}}), but this can
#'   be changed by the user with
#'   \code{\link[RcppParallel:setThreadOptions]{RcppParallel::setThreadOptions()}}.
#'
#'   An exception to the above is when this function is called within a
#'   \code{\link[foreach:foreach]{foreach}} parallel loop \strong{made by dtwclust}. If the parallel
#'   workers do not have the number of threads explicitly specified, this function will default to 1
#'   thread per worker.
#'
