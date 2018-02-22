ui <- tagList(
    shinyjs::useShinyjs(),
    navbarPage(
        title = "Time Series Clustering",
        # ==========================================================================================
        # Explore tab
        tabPanel(
            "Explore",
            sidebarLayout(
                # sidebar --------------------------------------------------------------------------
                sidebarPanel(
                    textInput(
                        "explore__ids",
                        label = h3("Integer IDs"),
                        value = "For example 1:10, c(1,5), etc."
                    ),
                    numericInput(
                        "explore__height",
                        label = h3("Plot height (px)"),
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
                # main -----------------------------------------------------------------------------
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
            # top ----------------------------------------------------------------------------------
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
            # clustering parameters ----------------------------------------------------------------
            fluidRow(
                column(
                    3L,
                    # shared -----------------------------------------------------------------------
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
                    # extra arguments --------------------------------------------------------------
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
                    # tsclust controls -------------------------------------------------------------
                    hr(),
                    h3("Control parameters"),
                    # partitional controls ---------------------------------------------------------
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
                    # hierarchical controls --------------------------------------------------------
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
                    # fuzzy controls ---------------------------------------------------------------
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
                    # tadpole controls -------------------------------------------------------------
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
