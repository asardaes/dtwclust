#' Note: any variable whose name starts with . and ends with _ was defined in the enclosing
#' environment of server when calling interactive_clustering().

server <- function(input, output, session) {
    # ==============================================================================================
    # reactive values
    result <- reactiveVal(NA)
    pair_ids <- reactiveVal(NA)
    best_window <- reactiveVal(NA_integer_)
    constraints <- reactiveVal(data.frame())
    window_flags <- reactiveVal(data.frame())
    # non-reactive values
    pair_tracker <- NA
    cluster_ids <- NA
    # functions
    score_fun <- function(obj_list, ...) {
        df_list <- lapply(obj_list, function(obj) {
            df <- as.data.frame(rbind(obj@cluster))
            colnames(df) <- paste0(".cl_", 1L:length(obj@cluster))
            df
        })
        dplyr::bind_rows(df_list)
    }
    enable_buttons <- function() {
        shinyjs::enable("cluster__must_link")
        shinyjs::enable("cluster__cannot_link")
        shinyjs::enable("cluster__dont_know")
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
            enable_buttons()
            pair_ids(new_pair)
        }
    }
    majority <- function(x) {
        ux <- sort(unique(x))
        ux[which.max(tabulate(match(x, ux)))]
    }
    if (is.null(complexity))
        complexity <- function(flags) {
            if (length(flags) <= 1L) return(0)
            sign_changes <- sum(abs(diff(flags)))
            max_consecutive_true <- rle(flags)
            if (any(max_consecutive_true$values))
                max_consecutive_true <- max(max_consecutive_true$lengths[max_consecutive_true$values])
            else
                max_consecutive_true <- 0
            # return
            cmp <- sign_changes / (length(flags) - 1L) / max_consecutive_true
            if (is.na(cmp)) cmp <- Inf
            cmp
        }
    feedback_handler <- function(pair, type) {
        cp <- switch(
            type,
            "must_link" = `==`,
            "cannot_link" = `!=`
        )
        df <- data.frame(
            series1 = pair[1L],
            series2 = pair[2L],
            link_type = type,
            stringsAsFactors = FALSE
        )
        flags <- logical(nrow(cluster_ids))
        i <- 1L
        while (i <= nrow(cluster_ids)) {
            # + 2L due to cluster_ids having config_id and window_size as first two columns
            flags[i] <- cp(cluster_ids[i, pair[1L] + 2L],
                           cluster_ids[i, pair[2L] + 2L])
            i <- i + 1L
        }
        df$complexity <- complexity(flags)
        df$best_window <- cluster_ids[which.max(flags), 2L]
        # update reactive values
        constraints(rbind(constraints(), df))
        window_flags({
            rbind(window_flags(), data.frame(
                constraint = paste0(type, "=", paste0(pair, collapse = ",")),
                window_size = cluster_ids$window_size,
                flag = as.integer(flags),
                stringsAsFactors = FALSE
            ))
        })
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
    observe({
        if (input$cluster__part_nrep > 1L) {
            shinyjs::enable("cluster__part_agg")
        }
        else {
            shinyjs::disable("cluster__part_agg")
        }
    })
    observe({
        if (!input$cluster__hier_method_custom && input$cluster__hier_method == "all") {
            shinyjs::enable("cluster__hier_agg")
        }
        else {
            shinyjs::disable("cluster__hier_agg")
        }
    })
    # ----------------------------------------------------------------------------------------------
    # cluster
    observeEvent(input$cluster__cluster, main, handler.quoted = TRUE)
    # main plot
    observeEvent(c(pair_ids(), input$cluster__plot_height), cluster_plot, handler.quoted = TRUE)
    # ----------------------------------------------------------------------------------------------
    # annotation feedback
    observe({
        cnst <- constraints()
        if (nrow(cnst) > 0L) {
            threshold <- input$cluster__complexity
            df <- dplyr::filter(cnst, complexity > 0 & complexity < threshold)
            trivial <- all(sapply(cnst$complexity, function(cx) {
                isTRUE(all.equal(cx, 0)) | is.infinite(cx)
            }))
            if (nrow(df) > 0L)
                best_window(majority(df$best_window))
            else if (trivial)
                best_window(min(cnst$best_window))
            else
                best_window(NA)
        }
    })
    observeEvent(input$cluster__must_link, {
        ids <- pair_ids()
        connected <- pair_tracker$link(ids[1L], ids[2L], 1L)
        feedback_handler(ids, "must_link")
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
        complete <- pair_tracker$link(ids[1L], ids[2L], 0L)
        feedback_handler(ids, "cannot_link")
        if (complete) {
            pair_ids(NULL)
            disable_buttons()
            if (complete)
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
        complete <- pair_tracker$link(ids[1L], ids[2L], -1L)
        if (complete) {
            pair_ids(NULL)
            disable_buttons()
            if (complete)
                shinyjs::alert(paste(
                    "No unlinked pairs left.",
                    "Based on your feedback,",
                    "nothing can be inferred."
                ))
        }
        else {
            get_new_pair()
        }
    })
    output$cluster__best_window <- renderText({
        paste("Best window size so far:", best_window())
    })
    # ==============================================================================================
    # Evaluate tab
    # ----------------------------------------------------------------------------------------------
    # summary
    output$evaluate__summary <- renderText({
        cnst <- constraints()
        out <- paste0("Suggested window: ", best_window(), "<br>\n",
                      "Annotations so far: ", nrow(cnst), "<br>\n")
        if (nrow(cnst) > 0L) {
            out <- paste0(
                out,
                "Number of 'must link': ", sum(cnst$link_type == "must_link"), "<br>\n",
                "Number of 'cannot link': ", sum(cnst$link_type == "cannot_link"), "<br>\n"
            )
        }
        out
    })
    # ----------------------------------------------------------------------------------------------
    # save
    observeEvent(input$evaluate__save, {
        res <- result()
        if (inherits(res, "list")) {
            res$ensembles <- cluster_ids
            res$constraints <- constraints()
            res$constraints_plot_df <- window_flags()
            res$best_window <- best_window()
            out_name <- input$evaluate__save_name
            if (nzchar(out_name)) {
                tryCatch({
                    assign(out_name, res, globalenv())
                    shinyjs::alert("Saved! Exit shiny app to update the global environment.")
                },
                error = function(e) {
                    shinyjs::alert(paste("Could not save:", e$message))
                })
            }
        }
    })
    # ----------------------------------------------------------------------------------------------
    # constraints table
    output$evaluate__constraints <- renderTable(constraints())
    # ----------------------------------------------------------------------------------------------
    # constraints plots
    observe(constraints_plot, quoted = TRUE)
}
