library(shiny)

#' Note: any variable whose name starts with . and ends with _ was defined in the enclosing
#' environment of server when calling interactive_clustering().

ui <- tagList(
    shinyjs::useShinyjs(),
    navbarPage(
        title = "Time Series Clustering",
        # ==========================================================================================
        # Explore tab
        tabPanel(
            "Explore",
            sidebarLayout(
                sidebarPanel(
                    textInput(
                        "explore__ids",
                        label = h3("Integer IDs"),
                        value = "For example 1:10, c(1,5), etc."
                    ),
                    numericInput(
                        "explore__height",
                        label = h3("Plot height (per facet, in pixels)"),
                        value = 600,
                        step = 10
                    ),
                    actionButton(
                        "explore__trigger_plot",
                        label = "Plot series"
                    ),
                    helpText(
                        "It may take a few seconds for the plot to be displayed."
                    )
                ),
                mainPanel(
                    plotOutput(outputId = "explore__plot")
                )
            )
        )
    )
)

server <- function(input, output) {
    # ==============================================================================================
    # Explore tab, plot
    observeEvent(input$explore__trigger_plot, {
        shinyjs::disable("explore__trigger_plot")
        height <- as.integer(input$explore__height) * NCOL(.series_[[1L]])
        output$explore__plot <- renderPlot(
            isolate({
                # SHINY-interactive-clustering.R
                explore__plot(input$explore__ids, .explore_df_, .series_)
            }),
            height = height
        )
        shinyjs::enable("explore__trigger_plot")
    })
}
