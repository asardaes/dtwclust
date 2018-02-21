#' Note: any variable whose name starts with . and ends with _ was defined in the enclosing
#' environment of server when calling interactive_clustering().

server <- function(input, output, session) {
    # ==============================================================================================
    # reactive values
    result <- reactiveVal(NA)
    pair_ids <- reactiveVal(NA)
    # non-reactive values
    pair_tracker <- NA
    score_fun <- function(obj_list, ...) {
        df_list <- lapply(obj_list, function(obj) {
            df <- as.data.frame(rbind(obj@cluster))
            colnames(df) <- paste0(".cl_", 1L:length(obj@cluster))
            df
        })
        dplyr::bind_rows(df_list)
    }
    disable_buttons <- function() {
        shinyjs::disable("cluster__must_link")
        shinyjs::disable("cluster__cannot_link")
        shinyjs::disable("cluster__dont_know")
    }
    get_new_pair <- function() {
        new_pair <- pair_tracker$get_unseen_pair()
        if (is.null(new_pair)) {
            pair_ids(NULL)
            disable_buttons()
            shinyjs::alert("No unlinked pairs left.")
        }
        else {
            pair_ids(new_pair)
        }
    }
    # ==============================================================================================
    # Explore tab
    # characteristics table
    output$explore__characteristics <- renderTable(characteristics_table, quoted = TRUE)
    # plot
    observeEvent(input$explore__trigger_plot, explore_plot, handler.quoted = TRUE)
    # ==============================================================================================
    # Cluster tab
    # ----------------------------------------------------------------------------------------------
    # distance options
    observe({
        if (input$cluster__clus_type == "t") {
            shinyjs::disable("cluster__dist_args")
        }
        else {
            shinyjs::enable("cluster__dist_args")
        }
    })
    # ----------------------------------------------------------------------------------------------
    # centroid options
    observe({
        if (input$cluster__cent_custom) {
            shinyjs::show("cluster__cent_func")
            shinyjs::hide("cluster__cent")
        }
        else {
            shinyjs::hide("cluster__cent_func")
            shinyjs::show("cluster__cent")
        }
    })
    observe({
        choices <- switch(
            input$cluster__clus_type,
            "p" = as.list(centroids_nonfuzzy),
            "h" = list("default"),
            "t" = list("default")
        )
        selected <- switch(
            input$cluster__clus_type,
            "p" = "pam",
            "h" = "default",
            "t" = "default"
        )
        updateSelectInput(
            session,
            "cluster__cent",
            label = "Centroid",
            choices = choices,
            selected = selected
        )
    })
    # ----------------------------------------------------------------------------------------------
    # control options
    observe({
        if (!input$cluster__cent_custom && input$cluster__cent == "pam") {
            shinyjs::enable("cluster__part_pam")
        }
        else {
            shinyjs::disable("cluster__part_pam")
        }
    })
    observe({
        if (input$cluster__hier_method_custom) {
            shinyjs::show("cluster__hier_method_func")
            shinyjs::hide("cluster__hier_method")
        }
        else {
            shinyjs::hide("cluster__hier_method_func")
            shinyjs::show("cluster__hier_method")
        }
    })
    # ----------------------------------------------------------------------------------------------
    # cluster
    observeEvent(input$cluster__cluster, main, handler.quoted = TRUE)
    # main plot
    observeEvent(c(pair_ids(), input$cluster__plot_height), cluster_plot, handler.quoted = TRUE)
    # ----------------------------------------------------------------------------------------------
    # annotation feedback
    observeEvent(input$cluster__must_link, {
        ids <- pair_ids()
        connected <- pair_tracker$link(ids[1L], ids[2L], 1L)
        if (connected) {
            pair_ids(NULL)
            disable_buttons()
            if (connected)
                shinyjs::alert(paste(
                    "No unlinked pairs left.",
                    "Based on your feedback,",
                    "all series should go in 1 cluster."
                ))
        }
        else {
            get_new_pair()
        }
    })
    observeEvent(input$cluster__cannot_link, {
        ids <- pair_ids()
        connected <- pair_tracker$link(ids[1L], ids[2L], 0L)
        if (connected) {
            pair_ids(NULL)
            disable_buttons()
            if (connected)
                shinyjs::alert(paste(
                    "No unlinked pairs left.",
                    "Based on your feedback,",
                    "each series should go in its own cluster."
                ))
        }
        else {
            get_new_pair()
        }
    })
    observeEvent(input$cluster__dont_know, {
        ids <- pair_ids()
        pair_tracker$link(ids[1L], ids[2L], -1L)
        get_new_pair()
    })
    # ==============================================================================================
    # Evaluate tab
    # ----------------------------------------------------------------------------------------------
    # summary
    output$evaluate__summary <- renderText({
        out <- ""
        res <- result()
        if (inherits(res, "list")) {
        }
        out
    })
    # ----------------------------------------------------------------------------------------------
    # save
    observeEvent(input$evaluate__save, {
        if (inherits(result(), "list")) {
            out_name <- input$evaluate__save_name
            if (nzchar(out_name)) {
                tryCatch({
                    assign(out_name, result(), globalenv())
                    shinyjs::alert("Saved! Exit shiny app to update the global environment.")
                },
                error = function(e) {
                    shinyjs::alert(paste("Could not save:", e$message))
                })
            }
        }
    })
}
