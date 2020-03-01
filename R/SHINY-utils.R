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
    dfm$L1 <- factor(dfm$L1, levels = names(series))
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
    if (!is.integer(ids)) ids <- as.integer(eval(parse(n = 1L, text = ids)))
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
