#' A shiny app for semi-supervised DTW-based clustering
#'
#' Display a shiny user interface that implements the approach in Dau et al. (2016).
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom shiny shinyApp
#' @importFrom shinyjs useShinyjs
#'
#' @param series Time series in the formats accepted by [compare_clusterings()].
#' @param ... More arguments for [shiny::runApp()].
#' @param complexity A function to calculate a constraint's complexity. See details in the Cluster
#'   section.
#'
#' @details
#'
#' The approach developed in Dau et al. (2016) argues that finding a good value of `window.size` for
#' the DTW distance is very important, and suggests how to find one by using user-provided feedback.
#' After clustering is done, a pair of series is presented at a time, and the user must annotate the
#' pair as:
#'
#' - Must link: the series should be in the same cluster.
#' - Cannot link: the series should *not* be in the same cluster.
#' - Skip: the choice is unclear.
#'
#' After each step, a good value of the window size is suggested by evaluating which clusterings
#' fulfill the constraint(s) so far, and how (see Dau et al. (2016) for more information), and
#' performing a majority vote using the window sizes inferred from each constraint. The (main)
#' procedure is thus interactive and can be abandoned at any point.
#'
#' @section Explore:
#'
#'   This part of the app is simply to see some basic characteristics of the provided series and
#'   plot some of them. The field for integer IDs expects a valid R expression that specifies which
#'   of the `series` should be plotted. Multivariate series are plotted with each variable in a
#'   different facet.
#'
#' @section Cluster:
#'
#'   This part of the app implements the main procedure by leveraging [compare_clusterings()]. The
#'   interface is similar to [interactive_clustering()], so it's worth checking its documentation
#'   too. Since [compare_clusterings()] supports parallelization with [foreach::foreach()], you can
#'   register a parallel backend before opening the shiny app, but you should pre-load the workers
#'   with the necessary packages and/or functions. See [parallel::clusterEvalQ()] and
#'   [parallel::clusterExport()], as well as the examples below.
#'
#'   The range of window sizes is specified with a slider, and represents the size as a percentage
#'   of the shortest series' length. The `step` parameter indicates how spaced apart should the
#'   sizes be (parameter `'by'` in [base::seq()]). A 0-size window should only be used if all series
#'   have the same length. If the series have different lengths, using small window sizes can be
#'   problematic if the length differences are very big, see the notes in [dtw_basic()].
#'
#'   A `window.size` should *not* be specified in the extra parameters, it will be replaced with the
#'   computed values based on the slider. Using [dba()] centroid is detected, and will use the same
#'   window sizes.
#'
#'   For partitional clusterings with many repetitions, and hierarchical clusterings with many
#'   linkage methods, the resulting partitions are aggregated by calling [clue::cl_medoid()] with
#'   the specified aggregation `method`.
#'
#'   By default, complexity of a constraint is calculated differently from what is suggested in Dau
#'   et al. (2016):
#'
#'   - Allocate a logical flag vector with length equal to the number of tested window sizes.
#'   - For each window size, set the corresponding flag to `TRUE` if the constraint given by the
#'     user is fulfilled.
#'   - Calculate complexity as: (number of sign changes in the vector) / (number of window sizes -
#'     1L) / (maximum number of *contiguous* `TRUE` flags).
#'
#'   You can provide your own function in the `complexity` parameter. It will receive the flag
#'   vector as only input, and a single number is expected as a result.
#'
#'   The complexity threshold can be specified in the app. Any constraint whose complexity is higher
#'   than the threshold will not be considered for the majority vote. Constraints with a complexity
#'   of 0 are also ignored. An infinite complexity means that the constraint is never fulfilled by
#'   any clustering.
#'
#' @section Evaluate:
#'
#'   This section provides numerical results for reference. The latest results can be saved in the
#'   global environment, which includes clustering results, constraints so far, and the suggested
#'   window size. Since this includes everything returned by [compare_clusterings()], you could also
#'   use [repeat_clustering()] afterwards.
#'
#'   The constraint plots depict if the constraints are fulfilled or not for the given window sizes,
#'   where 1 means it was fulfilled and 0 means it wasn't. An error about a zero-dimension viewport
#'   indicates the plot height is too small to fit the plots, so please increase the height.
#'
#' @note
#'
#' The optimization mentioned in section 3.4 of Dau et al. (2016) is also implemented here.
#'
#' Tracing is printed to the console.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Dau, H. A., Begum, N., & Keogh, E. (2016). Semi-supervision dramatically improves time series
#' clustering under dynamic time warping. In Proceedings of the 25th ACM International on Conference
#' on Information and Knowledge Management (pp. 999-1008). ACM.
#' \url{https://sites.google.com/site/dtwclustering/}
#'
#' @seealso
#'
#' [interactive_clustering()], [compare_clusterings()]
#'
#' @examples
#'
#' \dontrun{
#' require(doParallel)
#' workers <- makeCluster(detectCores())
#' clusterEvalQ(workers, {
#'     library(dtwclust)
#'     RcppParallel::setThreadOptions(1L)
#' })
#' registerDoParallel(workers)
#' ssdtwclust(reinterpolate(CharTraj[1L:20L], 150L))
#' }
#'
ssdtwclust <- function(series, ..., complexity = NULL) {
    series <- tslist(series)
    check_consistency(series, "vltslist")
    is_multivariate(series) # dimension consistency check
    if (!is.null(complexity) && !is.function(complexity))
        stop("A custom complexity should be a function.")

    file_path <- system.file("ssdtwclust", "app.R", package = "dtwclust")
    if (!nzchar(file_path)) stop("Shiny app not found")
    ui <- server <- NULL # avoid NOTE about undefined globals
    source(file_path, local = TRUE, chdir = TRUE)

    server_env <- environment(server) # will see all dtwclust functions
    server_env$.series_ <- series
    server_env$.explore_df_ <- explore__tidy_series(series)
    server_env$complexity <- complexity

    app <- shiny::shinyApp(ui, server)
    shiny::runApp(app, ...)
}
