# Explore tab, main plot
explore_plot <- quote({
    shinyjs::disable("explore__trigger_plot")
    height <- as.integer(input$explore__height)
    output$explore__plot <- renderPlot(
        isolate({
            # SHINY-utils.R
            explore__plot(input$explore__ids, .explore_df_, .series_)
        }),
        height = height
    )
    shinyjs::enable("explore__trigger_plot")
})

# Cluster tab, main plot
cluster_plot <- quote({
    if (inherits(result(), "list")) {
        tried <- tryCatch({
            ids <- pair_ids()
            if (is.null(ids)) {
                output$cluster__plot <- renderPlot({})
            }
            else {
                output$cluster__plot <- renderPlot(isolate({
                    explore__plot(ids, .explore_df_, .series_)
                }),
                height = as.integer(input$cluster__plot_height))
            }
        },
        error = function(e) {
            e
        })
        if (inherits(tried, "error")) {
            shinyjs::alert(tried$message)
        }
    }
    else {
        output$cluster__plot <- renderPlot({})
    }
})

# Evaluate tab, constraints plot
constraints_plot <- quote({
    df <- window_flags()
    cnst <- constraints()
    output$evaluate__plot <- renderPlot({
        if (nrow(df) == 0L) {
            NULL
        }
        else {
            df_labels <- data.frame(
                constraint = unique(df$constraint),
                complexity = paste("complexity =", round(cnst$complexity, 2L)),
                x = min(cluster_ids$window_size),
                y = 0.5,
                stringsAsFactors = FALSE
            )
            ggplot2::ggplot(df, ggplot2::aes(x = window_size, y = flag)) +
                ggplot2::geom_step(size = 2) +
                ggrepel::geom_label_repel(ggplot2::aes(x = x, y = y, label = complexity),
                                          data = df_labels,
                                          inherit.aes = FALSE,
                                          size = 10L) +
                ggplot2::facet_wrap(~constraint, ncol = 1L) +
                ggplot2::theme_bw(base_size = 20L)
        }
    },
    height = as.integer(input$evaluate__plot_height) * nrow(cnst))
})
