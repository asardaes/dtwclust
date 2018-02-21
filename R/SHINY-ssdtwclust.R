#' A shiny app for semi-supervised DTW-based clustering
#'
#' Display a shiny user interface that implements the approach in Dau et al. (2016).
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom shiny shinyApp
#' @importFrom shinyjs useShinyjs
#'
#' @param series Time series in the formats accepted by [compare_clusterings()].
#' @param ... More arguments for [shiny::runApp()].
#'
#' @details
#'
#' The approach suggested in Dau et al. (2016) argues that finding a good value of `window.size` for
#' the DTW distance is very important, and presents how to find one by using user-provided feedback.
#' After clustering is done, a pair of series is presented at a time, and the user must annotate the
#' pair as:
#'
#' - Must link: the series should be in the same cluster.
#' - Cannot link: the series should *not* be in the same cluster.
#' - Skip: the choice is unclear.
#'
#' After each step, a good value of the window size is suggested by evaluating which clusterings
#' fulfill the constraint(s) so far, and how (see Dau et al. (2016) for more information). The
#' (main) procedure is thus interactive and can be abandoned at any point.
#'
#' @section Explore:
#'
#'   This part of the app is simply to see some basic characteristics of the provided series and
#'   plot some of them. The field for integer IDs expects a valid R expression that specifies which
#'   of the `series` should be plotted. Multivariate series are plotted with each variable in a
#'   different facet.
#'
#' @section Cluster:
#'
#'   This part of the app implements the main procedure by leveraging [compare_clusterings()]. The
#'   interface is similar to [interactive_clustering()], so it is worth checking its documentation
#'   too. Since [compare_clusterings()] supports parallelization with [foreach::foreach()], you can
#'   register a parallel backend before opening the shiny app, but you should pre-load the workers
#'   with the necessary packages and/or functions. See [parallel::clusterEvalQ()] and
#'   [parallel::clusterExport()], as well as the examples below.
#'
#'   The range of window sizes is specified with the slider, and represents the size as a percentage
#'   of the shortest series' length. The `step` parameter indicates how spaced should the sizes be
#'   (parameter `'by'` in [base::seq()]). The `window.size` should *not* be specified in the extra
#'   parameters, it will be replaced with the computed values based on the slider. Using [dba()]
#'   centroid is detected, and will use the same window sizes.
#'
#'   For partitional clusterings with many repetitions and hierarchical clusterings with many
#'   linkage methods, the resulting partitions can be aggregated by calling [clue::cl_medoid()] with
#'   the specified `method`.
#'
#' @section Evaluate:
#'
#'   This section provides numerical results for reference.
#'
#' @note
#'
#' The optimization mentioned in section 3.4 of Dau et al. (2016) is also implemented here.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Dau, H. A., Begum, N., & Keogh, E. (2016). Semi-supervision dramatically improves time series
#' clustering under dynamic time warping. In Proceedings of the 25th ACM International on Conference
#' on Information and Knowledge Management (pp. 999-1008). ACM.
#' \url{https://sites.google.com/site/dtwclustering/}
#'
#' @examples
#'
#' \dontrun{
#' require(doParallel)
#' workers <- makeCluster(detectCores())
#' clusterEvalQ(workers, {
#'     library(dtwclust)
#'     RcppParallel::setThreadOptions(1L)
#' })
#' registerDoParallel(workers)
#' ssdtwclust(CharTrajMV)
#' }
#'
ssdtwclust <- function(series, ...) {
    series <- tslist(series)
    check_consistency(series, "vltslist")
    is_multivariate(series) # dimension consistency check

    file_path <- system.file("ssdtwclust", "app.R", package = "dtwclust")
    if (!nzchar(file_path)) stop("Shiny app not found")
    ui <- server <- NULL # avoid NOTE about undefined globals
    source(file_path, local = TRUE, chdir = TRUE)

    server_env <- environment(server) # will see all dtwclust functions
    server_env$.series_ <- series
    server_env$.explore_df_ <- explore__tidy_series(series)

    app <- shiny::shinyApp(ui, server)
    shiny::runApp(app, ...)
}
