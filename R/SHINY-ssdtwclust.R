#' A shiny app for semi-supervised DTW-based clustering
#'
#' Display a shiny user interface that implements the approach in .
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom shiny shinyApp
#' @importFrom shinyjs useShinyjs
#'
#' @param series Time series in the formats accepted by [compare_clusterings()].
#' @param ... More arguments for [shiny::runApp()].
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
#'   The `window.size` parameter should *not* be specified in the extra parameters, it will be
#'   replaced with the computed values based on the slider. Using [dba()] centroid is detected, and
#'   will use the same window sizes.
#'
#'   Since [compare_clusterings()] supports parallelization with [foreach::foreach()], you can
#'   register a parallel backend before opening the shiny app, but you should pre-load the workers
#'   with the necessary packages and/or functions. See [parallel::clusterEvalQ()] and
#'   [parallel::clusterExport()], as well as the examples below.
#'
#' @section Evaluate:
#'
#' @author Alexis Sarda-Espinosa
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
