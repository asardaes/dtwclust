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
    if (inherits(result(), "TSClusters")) {
        # shinyjs::disable("cluster__plot_type")
        # shinyjs::disable("cluster__plot_clus")
        # shinyjs::disable("cluster__plot_labels")
        # shinyjs::disable("cluster__plot_height")
        tried <- tryCatch({
            height <- as.integer(input$cluster__plot_height)
            type <- input$cluster__plot_type
            clus <- as.integer(parse_input(input$cluster__plot_clus, "c"))
            labels <- input$cluster__plot_labels
            labels <- if (nzchar(labels)) parse_input(labels) else NULL
            output$cluster__plot <- renderPlot({
                plot_bool <- if (type == "dendrogram") TRUE else FALSE
                seed <- result()@seed
                if (length(seed) == 1L) set.seed(seed)
                out <- plot(
                    result(),
                    clus = clus,
                    plot = plot_bool,
                    type = type,
                    labels = labels
                )
                if (inherits(out, "ggplot")) {
                    out <- out + facet_wrap(~cl, scales = "free_y", ncol = 2L)
                    out <- ggplot2::ggplot_build(out)$plot
                }
                out
            },
            height = height)
        },
        error = function(e) {
            e
        })
        # shinyjs::enable("cluster__plot_type")
        # shinyjs::enable("cluster__plot_clus")
        # shinyjs::enable("cluster__plot_labels")
        # shinyjs::enable("cluster__plot_height")
        if (inherits(tried, "error")) {
            shinyjs::alert(tried$message)
        }
    }
})
