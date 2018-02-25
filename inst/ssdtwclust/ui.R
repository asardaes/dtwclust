ui <- tagList(
    shinyjs::useShinyjs(),
    navbarPage(
        title = "Semi-Supervised DTW Clustering",
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
                )
            ),
            # clustering parameters ----------------------------------------------------------------
            fluidRow(
                column(
                    3L,
                    h2("Clustering parameters"),
                    # shared -----------------------------------------------------------------------
                    selectInput(
                        "cluster__clus_type",
                        label = "Type",
                        choices = list(
                            "Partitional" = "p",
                            "Hierarchical" = "h",
                            "TADPole" = "t"
                        )
                    ),
                    numericInput(
                        "cluster__k",
                        label = "Desired number of clusters",
                        value = 2,
                        min = 2
                    ),
                    fluidRow(
                        column(
                            9L,
                            sliderInput(
                                "cluster__windows",
                                label = "Range of window sizes (%)",
                                min = 0,
                                max = 100,
                                step = 1,
                                value = c(1, 20)
                            )
                        ),
                        column(
                            3L,
                            numericInput(
                                "cluster__windows_step",
                                "Step (%, internal)",
                                value = 1,
                                min = 1
                            )
                        )
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
                        label = "Distance parameters (except window.size)",
                        value = "",
                        width = "100%"
                    ),
                    textInput(
                        "cluster__cent_args",
                        label = "Centroid parameters (except window.size)",
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
                        fluidRow(
                            column(
                                6L,
                                numericInput(
                                    "cluster__part_iter",
                                    label = "Max. number of iterations",
                                    value = 20,
                                    min = 1
                                )
                            ),
                            column(
                                6L,
                                numericInput(
                                    "cluster__part_nrep",
                                    label = "Number of repetitions",
                                    value = 1,
                                    min = 1
                                )
                            )
                        ),
                        shinyjs::disabled(selectInput(
                            "cluster__part_agg",
                            label = "Aggregation method",
                            choices = list(
                                "euclidean",
                                "manhattan",
                                "comemberships",
                                "symdiff",
                                "Rand",
                                "GV1",
                                "BA/A",
                                "BA/D",
                                "BA/E",
                                "VI"
                            )
                        )),
                        checkboxInput(
                            "cluster__part_pam",
                            label = "pam.precompute",
                            value = TRUE
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
                                        "all",
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
                        selectInput(
                            "cluster__hier_agg",
                            label = "Aggregation method",
                            choices = list(
                                "euclidean",
                                "manhattan",
                                "comemberships",
                                "symdiff",
                                "Rand",
                                "GV1",
                                "BA/A",
                                "BA/D",
                                "BA/E",
                                "VI"
                            )
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
                                selectInput(
                                    "cluster__tadp_lb",
                                    label = "Lower bound",
                                    choices = list(
                                        "LB_Keogh" = "lbk",
                                        "LB_Improved" = "lbi"
                                    )
                                )
                            )
                        )
                    )
                ),
                # ----------------------------------------------------------------------------------
                # output plot
                column(
                    9L,
                    fluidRow(
                        column(
                            3L, offset = 4L,
                            fluidRow(
                                column(
                                    4L,
                                    shinyjs::disabled(actionButton(
                                        "cluster__must_link",
                                        label = "Must link"
                                    ))
                                ),
                                column(
                                    4L,
                                    shinyjs::disabled(actionButton(
                                        "cluster__cannot_link",
                                        label = "Cannot link"
                                    ))
                                ),
                                column(
                                    3L, offset = 1L,
                                    shinyjs::disabled(actionButton(
                                        "cluster__dont_know",
                                        label = "Skip"
                                    ))
                                )
                            )
                        ),
                        column(
                            2L, offset = 3L,
                            numericInput(
                                "cluster__plot_height",
                                label = "Plot height (px)",
                                value = 600,
                                step = 10
                            )
                        )
                    ),
                    fluidRow(
                        fluidRow(
                            column(
                                2L,
                                textOutput("cluster__best_window")
                            ),
                            column(
                                3L, offset = 2L,
                                sliderInput(
                                    "cluster__complexity",
                                    label = "Complexity threshold",
                                    min = 0,
                                    max = 1,
                                    value = 0.3,
                                    step = 0.01
                                )
                            )
                        ),
                        column(12L,
                               plotOutput(outputId = "cluster__plot")
                        )
                    )
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
                    helpText("The latest results",
                             "can be saved in the current R session's",
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
                            "Clustering results",
                            tableOutput("evaluate__raw")
                        ),
                        tabPanel(
                            "Constraints",
                            tableOutput("evaluate__constraints")
                        ),
                        tabPanel(
                            "Constraint plots",
                            numericInput(
                                "evaluate__plot_height",
                                label = "Plot height (per constraint, px)",
                                value = 200,
                                step = 10
                            ),
                            plotOutput("evaluate__plot")
                        )
                    )
                )
            )
        )
    )
)
