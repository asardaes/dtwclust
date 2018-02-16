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
interactive_clustering <- function(series, ...) {
    series <- tslist(series)
    check_consistency(series, "vltslist")
    is_multivariate(series) # dimension consistency check

    file_path <- system.file("interactive-clustering/app.R", package = "dtwclust")
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
    # add vertical lines to separate variables of multivariate series
    if (is_multivariate(series)) {
        gg <- gg + ggplot2::facet_wrap(~L2, ncol = 1L)
    }
    gg <- gg + ggplot2::theme_bw()
    # return
    gg
}
