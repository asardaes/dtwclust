#' A shiny app for interactive clustering
#'
#' Display a shiny user interface to do clustering based on the provided series.
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom shiny shinyApp
#' @importFrom shinyjs useShinyjs
#'
#' @param series Time series in the formats accepted by [tsclust()].
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
#'   This part of the app wraps [tsclust()], so you should be familiar with it. Some remarks:
#'
#'   - Specifying a custom centroid or hierarchical method expects the name of a function available
#'   in the R session (without quotes). Naturally, any required package should be loaded before
#'   calling `interactive_clustering`. For example, if you want to use [cluster::agnes()], you
#'   should load \pkg{cluster} beforehand.
#'   - A random seed of 0 means that it will be left as `NULL` when calling [tsclust()].
#'   - The input fields for Extra parameters (distance, centroid and ellipsis) expect a comma-
#'   separated sequence of key-value pairs. For example: `window.size = 10L`, `trace = TRUE`. You
#'   should be able to pass any variables available in the R session's global environment.
#'   - Regarding plot parameters:
#'     + The `Clusters` field is like the integer IDs from the Explore section.
#'     + The `Labels` field is passed to the plot method (see [TSClusters-methods]). You can specify
#'     several values like with the Extra parameters, e.g.: `nudge_x = 10`, `nudge_y = 1`. You can
#'     type an empty space to activate them with the defaults, and delete everything to hide them.
#'     Note that the location of the labels is random each time.
#'
#'   The plot area reacts to the plot parameters, but the actual clustering with [tsclust()] won't
#'   be executed until you click the `Cluster!` button. **The plot can take a couple of seconds to
#'   load!** Plotting multivariate series might generate warnings about missing values, they can be
#'   safely ignored.
#'
#'   Some of the control parameters are disabled when \pkg{dtwclust} detects them automatically.
#'
#'   The cross-distance matrix is cached so that it can be re-used when appropriate. The cached
#'   version is invalidated automatically when necessary.
#'
#' @section Evaluate:
#'
#'   This part of the app provides results of the current clustering. External CVIs can be
#'   calculated if the name of a variable with the ground truth is provided (see [cvi()]).
#'
#' @note
#'
#' Tracing is printed to the console.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @examples
#'
#' \dontrun{
#' interactive_clustering(CharTrajMV)
#' }
#'
interactive_clustering <- function(series, ...) {
    series <- tslist(series)
    check_consistency(series, "vltslist")
    is_multivariate(series) # dimension consistency check

    file_path <- system.file("interactive-clustering", "app.R", package = "dtwclust")
    if (!nzchar(file_path)) stop("Shiny app not found")
    ui <- server <- NULL # avoid NOTE about undefined globals
    source(file_path, local = TRUE, chdir = TRUE)

    server_env <- environment(server) # will see all dtwclust functions
    server_env$.series_ <- series
    server_env$.explore_df_ <- explore__tidy_series(series)

    app <- shiny::shinyApp(ui, server)
    shiny::runApp(app, ...)
}
