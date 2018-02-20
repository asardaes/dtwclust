# Explore tab, main plot
explore_plot <- quote({
    shinyjs::disable("explore__trigger_plot")
    height <- as.integer(input$explore__height) * NCOL(.series_[[1L]])
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
    if (inherits(result(), "list") && input$cluster__continue) {
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
