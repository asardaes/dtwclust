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
#'     in the R session (without quotes). Naturally, any required package should be loaded before
#'     calling `interactive_clustering`. For example, if you want to use [cluster::agnes()], you
#'     should load \pkg{cluster} beforehand.
#'   - A random seed of 0 means that it will be left as `NULL` when calling [tsclust()].
#'   - The input fields for Extra parameters (distance, centroid and ellipsis) expect a comma-
#'     separated sequence of key-value pairs. For example: `window.size = 10L`, `trace = TRUE`. You
#'     should be able to pass any variables available in the R session's global environment.
#'   - Regarding plot parameters:
#'     + The `Clusters` field is like the integer IDs from the Explore section.
#'     + The `Labels` field is passed to the plot method (see [TSClusters-methods]). You can specify
#'       several values like with the Extra parameters, e.g.: `nudge_x = 10`, `nudge_y = 1`. You can
#'       type an empty space to activate them with the defaults, and delete everything to hide them.
#'       Note that the location of the labels is random each time.
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
    if (!nzchar(file_path)) stop("Shiny application not found")
    ui <- server <- NULL # avoid NOTE about undefined globals
    source(file_path, local = TRUE, chdir = TRUE)

    server_env <- new.env() # will see all dtwclust functions
    server_env$.series_ <- series
    server_env$.explore_df_ <- explore__tidy_series(series)

    environment(server) <- server_env
    app <- shiny::shinyApp(ui, server)
    shiny::runApp(app, ...)
}

# ==================================================================================================
#' This helper will create the data frame used to plot in the Explore tab panel
#'
#' @keywords internal
#' @importFrom reshape2 melt
#'
#' @param series The list of time series
#'
#' @return A data frame
#'
explore__tidy_series <- function(series) {
    num_col <- NCOL(series[[1L]])
    num_rows <- sapply(series, NROW)
    if (num_col > 1L) {
        series <- lapply(series, function(this_series) {
            lapply(1L:num_col, function(j) { this_series[, j, drop = TRUE] })
        })
    }
    # transform to data frame
    dfm <- reshape2::melt(series)
    dfm$L1 <- factor(dfm$L1)
    if (num_col > 1L) dfm$L2 <- factor(dfm$L2)
    col_names <- colnames(dfm)
    col_names[col_names == "L1"] <- "Series"
    colnames(dfm) <- col_names
    # time indices
    dfm$t <- unlist(lapply(num_rows, function(len) {
        rep(seq(from = 1L, to = len), num_col)
    }))
    # return
    dfm
}

# ==================================================================================================
#' This helper will produce the plot in the Explore tab panel.
#'
#' @keywords internal
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#'
#' @param ids A character with the expression that should yield the integer ids.
#' @param df The data frame with the melted time series.
#' @param series The list of time series.
#'
explore__plot <- function(ids, df, series) {
    ids <- as.integer(eval(parse(n = 1L, text = ids)))
    df <- df[unclass(df$Series) %in% ids, , drop = FALSE]
    gg <- ggplot2::ggplot(df, ggplot2::aes_string(x = "t",
                                                  y = "value",
                                                  color = "Series",
                                                  group = "Series"))
    gg <- gg + ggplot2::geom_line()
    # add facets to separate variables of multivariate series
    if (is_multivariate(series)) {
        gg <- gg + ggplot2::facet_wrap(~L2, ncol = 1L)
    }
    gg <- gg + ggplot2::theme_bw()
    # return
    gg
}

# ==================================================================================================
#' This helper will parse comma-separated key-value pairs
#'
#' @keywords internal
#'
#' @param input_text A character.
#' @param into The name of the function that will contain the resulting values.
#'
parse_input <- function(input_text, into = "list") {
    input_text <- paste0(into, "(", input_text, ")")
    eval(parse(n = 1L, text = input_text))
}
