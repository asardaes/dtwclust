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
                    h4("Characteristics"),
                    tableOutput("explore__characteristics"),
                    h4("Plot"),
                    plotOutput(outputId = "explore__plot")
                )
            )
        ),
        # ==========================================================================================
        # Cluster tab
        tabPanel(
            "Cluster",
            # --------------------------------------------------------------------------------------
            # top
            fluidRow(
                column(
                    2L,
                    actionButton(
                        "cluster__cluster",
                        label = "Cluster!"
                    )
                ),
                column(
                    2L, offset = 8L,
                    h2("Plot parameters")
                )
            ),
            fluidRow(
                column(
                    3L,
                    h2("Clustering parameters")
                ),
                column(
                    2L, offset = 1L,
                    selectInput(
                        "cluster__plot_type",
                        label = "Type",
                        choices = list(
                            "sc" = 1L,
                            "series" = 2L,
                            "centroids" = 3L
                        )
                    )
                ),
                column(
                    2L,
                    textInput(
                        "cluster__plot_clus",
                        label = "Clusters",
                        value = "1:2"
                    )
                ),
                column(
                    2L,
                    textInput(
                        "cluster__plot_labels",
                        label = "Labels",
                        value = ""
                    )
                ),
                column(
                    2L,
                    numericInput(
                        "cluster__plot_height",
                        label = "Height (px)",
                        value = 600,
                        step = 10
                    )
                )
            ),
            # --------------------------------------------------------------------------------------
            # clustering parameters
            fluidRow(
                column(
                    3L,
                    selectInput(
                        "cluster__clus_type",
                        label = "Type",
                        choices = list(
                            "Partitional" = "p",
                            "Hierarchical" = "h",
                            "Fuzzy" = "f",
                            "TADPole" = "t"
                        )
                    ),
                    numericInput(
                        "cluster__k",
                        label = "Desired number of clusters",
                        value = 2,
                        min = 2
                    ),
                    selectInput(
                        "cluster__dist",
                        label = "Distance measure",
                        choices = with(new.env(), {
                            dist <- summary(proxy::pr_DB)
                            unname(lapply(dist$names[dist$distance], function(d) {
                                tolower(d[1L])
                            }))
                        }),
                        selected = "dtw_basic"
                    ),
                    fluidRow(
                        column(
                            9L,
                            selectInput(
                                "cluster__cent",
                                label = "Centroid",
                                choices = as.list(centroids_included)
                            ),
                            shinyjs::hidden(textInput(
                                "cluster__cent_func",
                                label = "Centroid",
                                value = "Function name"
                            ))
                        ),
                        column(
                            3L,
                            h5("", style = "padding:5px"),
                            checkboxInput(
                                "cluster__cent_custom",
                                label = "Custom",
                                value = FALSE
                            )
                        )
                    ),
                    numericInput(
                        "cluster__seed",
                        label = "Random seed",
                        value = 0,
                        min = 0
                    ),
                    checkboxInput(
                        "cluster__trace",
                        label = "Trace",
                        value = TRUE
                    ),
                    # ------------------------------------------------------------------------------
                    # extra arguments
                    hr(),
                    h3("Extra parameters"),
                    textInput(
                        "cluster__dist_args",
                        label = "Distance parameters",
                        value = "",
                        width = "100%"
                    ),
                    textInput(
                        "cluster__cent_args",
                        label = "Centroid parameters",
                        value = "",
                        width = "100%"
                    ),
                    textInput(
                        "cluster__dots",
                        label = "Ellipsis (...)",
                        value = "",
                        width = "100%"
                    ),
                    # ------------------------------------------------------------------------------
                    # tsclust controls
                    hr(),
                    h3("Control parameters"),
                    # partitional controls --------
                    conditionalPanel(
                        "input.cluster__clus_type == 'p'",
                        numericInput(
                            "cluster__part_iter",
                            label = "Max. number of iterations",
                            value = 20,
                            min = 1
                        ),
                        fluidRow(
                            column(
                                6L,
                                checkboxInput(
                                    "cluster__part_pam",
                                    label = "pam.precompute",
                                    value = TRUE
                                )
                            ),
                            column(
                                6L,
                                checkboxInput(
                                    "cluster__part_symmetric",
                                    label = "symmetric",
                                    value = FALSE
                                )
                            )
                        )
                    ),
                    # hierarchical controls --------
                    conditionalPanel(
                        "input.cluster__clus_type == 'h'",
                        fluidRow(
                            column(
                                9L,
                                selectInput(
                                    "cluster__hier_method",
                                    label = "method",
                                    choices = list(
                                        "ward.D",
                                        "ward.D2",
                                        "single",
                                        "complete",
                                        "average",
                                        "mcquitty",
                                        "median",
                                        "centroid"
                                    )
                                ),
                                shinyjs::hidden(textInput(
                                    "cluster__hier_method_func",
                                    label = "method",
                                    value = "Function name"
                                ))
                            ),
                            column(
                                3L,
                                h5("", style = "padding:5px"),
                                checkboxInput(
                                    "cluster__hier_method_custom",
                                    label = "Custom",
                                    value = FALSE
                                )
                            )
                        ),
                        checkboxInput(
                            "cluster__hier_symmetric",
                            label = "symmetric",
                            value = FALSE
                        )
                    ),
                    # fuzzy controls --------
                    conditionalPanel(
                        "input.cluster__clus_type == 'f'",
                        numericInput(
                            "cluster__fuzz_iter",
                            label = "Max. number of iterations",
                            value = 100,
                            min = 1
                        ),
                        fluidRow(
                            column(
                                6L,
                                numericInput(
                                    "cluster__fuzz_m",
                                    label = "Fuzziness",
                                    value = 2,
                                    min = 1
                                )
                            ),
                            column(
                                6L,
                                numericInput(
                                    "cluster__fuzz_delta",
                                    label = "Tolerance (delta)",
                                    value = 0.001,
                                    min = 0
                                )
                            )
                        ),
                        checkboxInput(
                            "cluster__fuzz_symmetric",
                            label = "symmetric",
                            value = FALSE
                        )
                    ),
                    # tadpole controls --------
                    conditionalPanel(
                        "input.cluster__clus_type == 't'",
                        fluidRow(
                            column(
                                6L,
                                numericInput(
                                    "cluster__tadp_dc",
                                    label = "Cutoff distance",
                                    value = 1,
                                    min = 0
                                )
                            ),
                            column(
                                6L,
                                numericInput(
                                    "cluster__tadp_window",
                                    label = "Window size",
                                    value = 10,
                                    min = 1
                                )
                            )
                        ),
                        selectInput(
                            "cluster__tadp_lb",
                            label = "Lower bound",
                            choices = list(
                                "LB_Keogh" = "lbk",
                                "LB_Improved" = "lbi"
                            )
                        )
                    )
                ),
                # ----------------------------------------------------------------------------------
                # output plot
                column(
                    9L,
                    plotOutput(outputId = "cluster__plot")
                )
            )
        ),
        # ==========================================================================================
        # Evaluate tab
        tabPanel(
            "Evaluate",
            sidebarLayout(
                sidebarPanel(
                    h3("Summary"),
                    htmlOutput("evaluate__summary"),
                    helpText("The latest TSClusters object can be saved in the current R session's",
                             "global environment by specifying the desired variable name",
                             "and clicking the save button."),
                    textInput(
                        "evaluate__save_name",
                        label = NULL
                    ),
                    actionButton(
                        "evaluate__save",
                        label = "Save"
                    )
                ),
                mainPanel(
                    tabsetPanel(
                        tabPanel(
                            "Clusters",
                            tableOutput("evaluate__clusinfo"),
                            tableOutput("evaluate__cl"),
                            tableOutput("evaluate__fcluster")
                        ),
                        tabPanel(
                            "CVIs",
                            fluidRow(
                                column(
                                    4L,
                                    textInput(
                                        "evaluate__ground_truth",
                                        label = "Ground truth variable"
                                    ),
                                    shinyjs::hidden(checkboxGroupInput(
                                        "evaluate__ext_cvis",
                                        label = "External Cluster Validity Indices",
                                        choices = list("external"),
                                        width = "100%"
                                    )),
                                    checkboxGroupInput(
                                        "evaluate__int_cvis",
                                        label = "Internal Cluster Validity Indices",
                                        choices = list("internal"),
                                        width = "100%"
                                    )
                                ),
                                column(
                                    8L,
                                    tableOutput("evaluate__cvis")
                                )
                            )
                        ),
                        tabPanel(
                            "Cross-distance matrix",
                            tableOutput("evaluate__distmat")
                        )
                    )
                )
            )
        )
    )
)

server <- function(input, output, session) {
    # ==============================================================================================
    # reactive values
    result <- reactiveVal(NA)
    distmat <- reactiveVal()
    # ==============================================================================================
    # Explore tab
    # ----------------------------------------------------------------------------------------------
    # characteristics table
    output$explore__characteristics <- renderTable({
        data.frame(
            "Amount of series" = length(.series_),
            "Multivariate" = is_multivariate(.series_),
            "Same lengths" = all(diff(lengths(.series_)) == 0L),
            "Centered" = all(sapply(.series_, function(series) {
                if (is.null(dim(series)))
                    mu <- mean(series)
                else
                    mu <- colMeans(series)
                all(sapply(mu, function(m) { isTRUE(all.equal(0, m)) }))
            })),
            "Scaled" = all(sapply(.series_, function(series) {
                if (is.null(dim(series)))
                    sigma <- sd(series)
                else
                    sigma <- apply(series, 2L, sd)
                all(sapply(sigma, function(s) { isTRUE(all.equal(1, s)) }))
            })),
            check.names = FALSE
        )
    })
    # ----------------------------------------------------------------------------------------------
    # plot(s)
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
    observeEvent(input$cluster__cluster, {
        shinyjs::disable("cluster__cluster")
        this_result <- tryCatch({
            type <- input$cluster__clus_type
            k <- as.integer(input$cluster__k)
            distance <- input$cluster__dist
            seed <- as.integer(input$cluster__seed)
            if (seed == 0L) seed <- NULL
            trace <- input$cluster__trace
            if (input$cluster__cent_custom)
                centroid <- match.fun(input$cluster__cent_func)
            else
                centroid <- input$cluster__cent
            # controls
            control <- switch(
                input$cluster__clus_type,
                "p" = {
                    partitional_control(
                        iter.max = as.integer(input$cluster__part_iter),
                        pam.precompute = input$cluster__part_pam,
                        symmetric = input$cluster__part_symmetric,
                        distmat = distmat()
                    )
                },
                "h" = {
                    if (input$cluster__hier_method_custom)
                        method <- match.fun(input$cluster__hier_method_func)
                    else
                        method <- input$cluster__hier_method

                    hierarchical_control(
                        method = method,
                        symmetric = input$cluster__hier_symmetric,
                        distmat = distmat()
                    )
                },
                "f" = {
                    fuzzy_control(
                        iter.max = as.integer(input$cluster__fuzz_iter),
                        fuzziness = input$cluster__fuzz_m,
                        delta = input$cluster__fuzz_delta,
                        symmetric = input$cluster__fuzz_symmetric,
                        distmat = distmat()
                    )
                },
                "t" = {
                    tadpole_control(
                        dc = input$cluster__tadp_dc,
                        window.size = as.integer(input$cluster__tadp_window),
                        lb = input$cluster__tadp_lb
                    )
                }
            )
            # dist args
            dist_args <- input$cluster__dist_args
            dist_args <- if (nzchar(dist_args)) parse_input(dist_args) else list()
            # cent args
            cent_args <- input$cluster__cent_args
            cent_args <- if (nzchar(cent_args)) parse_input(cent_args) else list()
            # dots
            dots <- input$cluster__dots
            dots <- if (nzchar(dots)) parse_input(dots) else list()
            # return
            args <- enlist(
                series = .series_,
                type = type,
                k = k,
                distance = distance,
                centroid = centroid,
                control = control,
                seed = seed,
                trace = trace,
                error.check = FALSE,
                args = tsclust_args(dist = dist_args, cent = cent_args),
                dots = dots
            )
            if (is.character(centroid) && centroid == "default")
                args$centroid <- NULL
            if (type == "t")
                args$distance <- NULL
            do.call(tsclust, args, TRUE)
        },
        error = function(e) {
            e
        })
        shinyjs::enable("cluster__cluster")
        if (inherits(this_result, "error")) {
            shinyjs::alert(this_result$message)
        }
        else {
            result(this_result)
            distmat(this_result@distmat)
        }
    })
    # ----------------------------------------------------------------------------------------------
    # main plot
    observe({
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
                    assign(out_name, result(), globalenv())
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
                output$evaluate__cvis <- renderTable({
                    as.data.frame(rbind(cvi(result(), b = b, type = c(internal, external))))
                })
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
                    out@fcluster,
                    rownames = TRUE,
                    caption = "Fuzzy memberships",
                    caption.placement = "top"
                )
            }
            else {
                output$evaluate__fcluster <- renderTable(data.frame())

                output$evaluate__clusinfo <- renderTable({
                    df <- t(out@clusinfo)
                    rownames(df) <- c("Cluster size", "Average distance")
                    colnames(df) <- paste("Cluster", 1L:out@k)
                    df
                },
                rownames = TRUE,
                caption = "Cluster sizes with average intra-cluster distance.",
                caption.placement = "top"
                )
                output$evaluate__cl <- renderTable({
                    df <- cbind(out@cluster, out@cldist)
                    colnames(df) <- c("Cluster ID", "Distance to centroid")
                    rownames(df) <- names(.series_)
                    df
                },
                rownames = TRUE,
                caption = "Cluster indices.",
                caption.placement = "top")
            }
        }
        if (is.null(dm))
            output$evaluate__distmat <- renderTable(data.frame())
        else
            output$evaluate__distmat <- renderTable(as.data.frame(dm),
                                                    spacing = "xs",
                                                    rownames = TRUE)
    })
}
