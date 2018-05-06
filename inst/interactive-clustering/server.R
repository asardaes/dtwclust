#' Note: any variable whose name starts with . and ends with _ was defined in the enclosing
#' environment of server when calling interactive_clustering().

server <- function(input, output, session) {
    # ==============================================================================================
    # reactive values
    result <- reactiveVal(NA)
    distmat <- reactiveVal()
    # ==============================================================================================
    # Explore tab
    # characteristics table
    output$explore__characteristics <- renderTable(characteristics_table, quoted = TRUE)
    # plot
    observeEvent(input$explore__trigger_plot, explore_plot, handler.quoted = TRUE)
    # ==============================================================================================
    # Cluster tab
    # ----------------------------------------------------------------------------------------------
    # invalidate cached distmat
    observeEvent(c(
        input$cluster__dist,
        input$cluster__dist_args,
        input$cluster__part_pam,
        input$cluster__part_symmetric,
        input$cluster__hier_symmetric,
        input$cluster__fuzz_symmetric
    ),
    {
        distmat(NULL)
    })
    observe({
        if (input$cluster__cent_custom) {
            distmat(NULL)
        }
        else if (!(input$cluster__cent %in% c("pam", "fcmdd", "default"))) {
            distmat(NULL)
        }
    })
    # ----------------------------------------------------------------------------------------------
    # plot type selection options
    observe({
        if (input$cluster__clus_type == "h")
            choices <- list("sc", "series", "centroids", "dendrogram")
        else
            choices <- list("sc", "series", "centroids")

        updateSelectInput(
            session,
            "cluster__plot_type",
            label = "Type",
            choices = choices
        )
    })
    observe({
        if (input$cluster__plot_type == "dendrogram") {
            shinyjs::disable("cluster__plot_clus")
            shinyjs::disable("cluster__plot_labels")
        }
        else {
            shinyjs::enable("cluster__plot_clus")
            shinyjs::enable("cluster__plot_labels")
        }
    })
    # ----------------------------------------------------------------------------------------------
    # distance options
    observe({
        if (input$cluster__clus_type == "t") {
            shinyjs::disable("cluster__dist")
            shinyjs::disable("cluster__dist_args")
        }
        else {
            shinyjs::enable("cluster__dist")
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
            "f" = as.list(centroids_fuzzy),
            "t" = list("default")
        )
        selected <- switch(
            input$cluster__clus_type,
            "p" = "pam",
            "h" = "default",
            "f" = "fcmdd",
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
        if (input$cluster__dist %in% distances_included) {
            shinyjs::disable("cluster__part_symmetric")
            shinyjs::disable("cluster__hier_symmetric")
            shinyjs::disable("cluster__fuzz_symmetric")
        }
        else {
            shinyjs::enable("cluster__part_symmetric")
            shinyjs::enable("cluster__hier_symmetric")
            shinyjs::enable("cluster__fuzz_symmetric")
        }
    })
    # ----------------------------------------------------------------------------------------------
    # cluster
    observeEvent(input$cluster__cluster, main, handler.quoted = TRUE)
    # main plot
    observe(cluster_plot, quoted = TRUE)
    # ==============================================================================================
    # Evaluate tab
    # ----------------------------------------------------------------------------------------------
    # summary
    output$evaluate__summary <- renderText({
        out <- ""
        tsc <- result()
        if (inherits(tsc, "TSClusters")) {
            dist <- tsc@distance
            isolate({ if (dist == "unknown") dist <- input$cluster__dist })
            out <- paste0(
                tools::toTitleCase(tsc@type),
                " clustering with ",
                tolower(dist),
                " distance and ",
                tolower(tsc@centroid),
                " centroids.<br>\n"
            )
            out <- paste0(
                out,
                switch(
                    tsc@type,
                    "partitional" =, "fuzzy" = {
                        converged <- if (tsc@converged) "(converged)" else "(did not converge)"
                        paste0(
                            tsc@k,
                            " clusters were created. It took ",
                            tsc@iter,
                            " iterations ",
                            converged,
                            ".<br>\n"
                        )
                    },
                    "hierarchical" = {
                        paste0(
                            "The ",
                            tsc@method,
                            " linkage was used, and ",
                            tsc@k,
                            " clusters were created with the cutree() function.<br>\n"
                        )
                    },
                    "tadpole" = {
                        lb <- switch(tsc@control$lb,
                                     "lbk" = "Keogh's",
                                     "lbi" = "Lemire's improved")
                        paste0(
                            tsc@k,
                            " clusters were created, using a cutoff distance of ",
                            tsc@control$dc,
                            " and a window size of size ",
                            tsc@control$window.size,
                            ". ",
                            lb,
                            " lower bound was used.<br>\n"
                        )
                    }
                )
            )
            out <- paste0(out,
                          "The measured execution time was ",
                          tsc@proctime[["elapsed"]],
                          " seconds.<br>\n")
        }
        out
    })
    # ----------------------------------------------------------------------------------------------
    # save
    observeEvent(input$evaluate__save, {
        if (inherits(result(), "TSClusters")) {
            out_name <- input$evaluate__save_name
            if (nzchar(out_name)) {
                tryCatch({
                    assign(out_name, result(), .GlobalEnv)
                    shinyjs::alert("Saved! Exit shiny app to update the global environment.")
                },
                error = function(e) {
                    shinyjs::alert(paste("Could not save:", e$message))
                })
            }
        }
    })
    # ----------------------------------------------------------------------------------------------
    # cvi
    # selection options
    observe({
        res <- result()
        if (inherits(res, "TSClusters")){
            if (inherits(res, "FuzzyTSClusters")) {
                if (nzchar(input$evaluate__ground_truth))
                    shinyjs::show("evaluate__ext_cvis")
                else
                    shinyjs::hide("evaluate__ext_cvis")
                # external
                choices <- list(
                    "Soft rand (max)" = "RI",
                    "Soft adjusted rand (max)" = "ARI",
                    "Soft variation of information (min)" = "VI",
                    "Soft normalized mutual information (max)" = "NMIM"
                )
                updateCheckboxGroupInput(
                    session,
                    "evaluate__ext_cvis",
                    label = "External Cluster Validity Indices",
                    choices = choices,
                    selected = isolate(input$evaluate__ext_cvis)
                )
                # internal
                choices <- list(
                    "MPC (max)" = "MPC",
                    "K (min)" = "K",
                    "T (min)" = "T",
                    "SC (max)" = "SC",
                    "PBMF (max)" = "PBMF"
                )
                updateCheckboxGroupInput(
                    session,
                    "evaluate__int_cvis",
                    label = "Internal Cluster Validity Indices",
                    choices = choices,
                    selected = isolate(input$evaluate__int_cvis)
                )
            }
            else {
                if (nzchar(input$evaluate__ground_truth))
                    shinyjs::show("evaluate__ext_cvis")
                else
                    shinyjs::hide("evaluate__ext_cvis")
                # external
                choices <- list(
                    "Rand (max)" = "RI",
                    "Adjusted rand (max)" = "ARI",
                    "Jaccard (max)" = "J",
                    "Fowlkes-Mallows (max)" = "FM",
                    "Variation of information (min)" = "VI"
                )
                updateCheckboxGroupInput(
                    session,
                    "evaluate__ext_cvis",
                    label = "External Cluster Validity Indices",
                    choices = choices,
                    selected = isolate(input$evaluate__ext_cvis)
                )
                # internal
                choices <- list(
                    "Silhouette (max)" = "Sil",
                    "Dunn (max)" = "D",
                    "COP (min)" = "COP",
                    "Davies-Bouldin (min)" = "DB",
                    "Modified Davies-Bouldin DB* (min)" = "DBstar",
                    "Calinski-Harabasz (max)" = "CH",
                    "Score function (max)" = "SF"
                )
                updateCheckboxGroupInput(
                    session,
                    "evaluate__int_cvis",
                    label = "Internal Cluster Validity Indices",
                    choices = choices,
                    selected = isolate(input$evaluate__int_cvis)
                )
            }
        }
    })
    # calcluation of cvi
    observe({
        if (inherits(result(), "TSClusters")) {
            internal <- input$evaluate__int_cvis
            truth <- input$evaluate__ground_truth
            tryCatch({
                if (nzchar(truth)) {
                    external <- input$evaluate__ext_cvis
                    b <- eval(parse(n = 1L, text = truth))
                }
                else {
                    external <- ""
                    b <- NULL
                }
                output$evaluate__cvis <- renderTable(cvis_table, quoted = TRUE)
            },
            error = function(e) {
                shinyjs::alert(e$message)
            })
        }
    })
    # ----------------------------------------------------------------------------------------------
    # clusters and distmat
    observe({
        out <- result()
        dm <- NULL
        if (inherits(out, "TSClusters")) {
            dm <- out@distmat

            if (inherits(out, "FuzzyTSClusters")) {
                output$evaluate__clusinfo <- renderTable(data.frame())
                output$evaluate__cl <- renderTable(data.frame())
                output$evaluate__fcluster <- renderTable(
                    fcluster_table,
                    quoted = TRUE,
                    rownames = TRUE,
                    caption = "Fuzzy memberships",
                    caption.placement = "top"
                )
            }
            else {
                output$evaluate__fcluster <- renderTable(data.frame())
                output$evaluate__clusinfo <- renderTable(
                    clusinfo_table,
                    quoted = TRUE,
                    rownames = TRUE,
                    caption = "Cluster sizes with average intra-cluster distance.",
                    caption.placement = "top"
                )
                output$evaluate__cl <- renderTable(
                    cl_table,
                    quoted = TRUE,
                    rownames = TRUE,
                    caption = "Cluster indices.",
                    caption.placement = "top"
                )
            }
        }
        if (is.null(dm))
            output$evaluate__distmat <- renderTable(data.frame())
        else
            output$evaluate__distmat <- renderTable(
                distmat_table,
                quoted = TRUE,
                spacing = "xs",
                rownames = TRUE
            )
    })
}
