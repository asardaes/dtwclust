#' @section Proxy version:
#'
#'   It also includes symmetric optimizations to calculate only half a distance matrix when
#'   appropriate---only one list of series should be provided in \code{x}. If you want to avoid this
#'   optimization, call \code{\link[proxy]{dist}} by giving the same list of series in both \code{x}
#'   and \code{y}.
#'
