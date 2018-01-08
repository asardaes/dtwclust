#' @section Parallel Computing:
#'
#'   Please note that running tasks in parallel does \strong{not} guarantee faster computations. The
#'   overhead introduced is sometimes too large, and it's better to run tasks sequentially.
#'
#'   This function uses the \code{\link[RcppParallel:RcppParallel-package]{RcppParallel}} package
#'   for parallelization. It uses all available threads by default, but this can be changed by the
#'   user with \code{\link[RcppParallel:setThreadOptions]{RcppParallel::setThreadOptions()}}. Note
#'   that the package is loaded but not attached by \pkg{dtwclust}.
#'
