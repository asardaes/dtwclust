#' Time series clustering along with optimizations for the Dynamic Time Warping distance
#'
#' Time series clustering with a wide variety of strategies and a series of optimizations specific
#' to the Dynamic Time Warping (DTW) distance and its corresponding lower bounds (LBs).
#'
#' @docType package
#' @name dtwclust-package
#'
#' @details
#'
#' Many of the algorithms implemented in this package are specifically tailored to DTW, hence its
#' name. However, the main clustering function is flexible so that one can test many different
#' clustering approaches, using either the time series directly, or by applying suitable
#' transformations and then clustering in the resulting space. Other implementations included in the
#' package provide some alternatives to DTW.
#'
#' DTW is a dynamic programming algorithm that tries to find the optimum warping path between two
#' series. Over the years, several variations have appeared in order to make the procedure faster or
#' more efficient. Please refer to the included references for more information, especially Giorgino
#' (2009), which is a good practical introduction.
#'
#' Most optimizations require equal dimensionality, which means time series should have equal
#' length. DTW itself does not require this, but it is relatively expensive to compute. Other
#' distance definitions may be used, or series could be reinterpolated to a matching length
#' (Ratanamahatana and Keogh 2004).
#'
#' The main clustering function and entry point for this package is [tsclust()], with a convenience
#' wrapper for multiple tests in [compare_clusterings()], and a shiny app in
#' [interactive_clustering()]. There is another less-general-purpose shiny app in [ssdtwclust()].
#'
#' Please note the random number generator is set to L'Ecuyer-CMRG when \pkg{dtwclust} is attached
#' in an attempt to preserve reproducibility. You are free to change this afterwards if you wish
#' (see [base::RNGkind()]), but \pkg{dtwclust} will always use L'Ecuyer-CMRG internally.
#'
#' For more information, please read the included package vignettes, which can be accessed by typing
#' `browseVignettes("dtwclust")`.
#'
#' @note
#'
#' This software package was developed independently of any organization or institution that is or
#' has been associated with the author.
#'
#' This package can be used without attaching it with [base::library()] with some caveats:
#'
#' - The \pkg{methods} [package][methods::methods-package] must be attached. `R` usually does this
#'   automatically, but [utils::Rscript()] does not.
#' - If you want to use the \pkg{proxy} version of [dtw::dtw()] (e.g. for clustering), you have to
#'   attach the \pkg{dtw} package manually.
#'
#' Be careful with reproducibility, `R`'s random number generator is only changed session-wide if
#' \pkg{dtwclust} is attached.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package's vignette's references.
#'
#' @seealso
#'
#' [tsclust()], [compare_clusterings()], [interactive_clustering()], [ssdtwclust()], [dtw_basic()],
#' [proxy::dist()].
#'
#' @useDynLib dtwclust, .registration = TRUE
#' @import foreach
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
NULL

# PREFUN for some of my proxy distances so that they support 'pairwise' directly
proxy_prefun <- function(x, y, pairwise, params, reg_entry) {
    params$pairwise <- pairwise
    list(x = x, y = y, pairwise = pairwise, p = params, reg_entry = reg_entry)
}

#' @importFrom utils packageVersion
#'
.onLoad <- function(lib, pkg) {
    # Register DTW2
    if (!check_consistency("DTW2", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw2_proxy, names=c("DTW2", "dtw2"),
                               loop = TRUE, type = "metric", distance = TRUE,
                               description = "DTW with L2 norm",
                               PACKAGE = "dtwclust")
    # Register DTW_BASIC
    if (!check_consistency("DTW_BASIC", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw_basic_proxy, names=c("DTW_BASIC", "dtw_basic"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Basic and maybe faster DTW distance",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register LB_Keogh
    if (!check_consistency("LB_Keogh", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_keogh_proxy, names=c("LBK", "LB_Keogh", "lbk"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Keogh's DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register LB_Improved
    if (!check_consistency("LB_Improved", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_improved_proxy, names=c("LBI", "LB_Improved", "lbi"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Lemire's improved DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register SBD
    if (!check_consistency("SBD", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = sbd_proxy, names=c("SBD", "sbd"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Paparrizos and Gravanos' shape-based distance for time series",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 2 - d })
    # Register DTW_LB
    if (!check_consistency("DTW_LB", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw_lb, names=c("DTW_LB", "dtw_lb"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "DTW distance aided with Lemire's lower bound",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register GAK
    if (!check_consistency("GAK", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = gak_proxy, names=c("GAK", "gak"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Fast (triangular) global alignment kernel distance",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 1 - d })
    # Register uGAK
    if (!check_consistency("uGAK", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = gak_simil, names=c("uGAK", "ugak"),
                               loop = FALSE, type = "metric", distance = FALSE,
                               description = "Fast (triangular) global alignment kernel similarity",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register soft-DTW
    if (!check_consistency("sdtw", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = sdtw_proxy, names=c("sdtw", "SDTW", "soft-DTW"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Soft-DTW",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    # avoids default message if no backend exists
    if (is.null(foreach::getDoParName())) foreach::registerDoSEQ()
}

.onAttach <- function(lib, pkg) {
    RNGkind(dtwclust_rngkind)

    packageStartupMessage("dtwclust:\n",
                          "Setting random number generator to L'Ecuyer-CMRG (see RNGkind()).\n",
                          'To read the included vignettes type: browseVignettes("dtwclust").\n',
                          'See news(package = "dtwclust") after package updates.')

    if (grepl("\\.9000$", utils::packageVersion("dtwclust")))
        packageStartupMessage("This is a developer version of dtwclust:\n",
                              "Using devtools::test() is currently broken, see tests/testthat.R")
}

.onUnload <- function(libpath) {
    # Unegister distances
    if (check_consistency("DTW2", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("DTW2")
    if (check_consistency("DTW_BASIC", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("DTW_BASIC")
    if (check_consistency("LB_Keogh", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("LB_Keogh")
    if (check_consistency("LB_Improved", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("LB_Improved")
    if (check_consistency("SBD", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("SBD")
    if (check_consistency("DTW_LB", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("DTW_LB")
    if (check_consistency("GAK", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("GAK")
    if (check_consistency("uGAK", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("uGAK")
    if (check_consistency("sdtw", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("sdtw")
    library.dynam.unload("dtwclust", libpath)
}

release_questions <- function() {
    c(
        "Changed .Rbuildignore to exclude test rds files?",
        "Built the binary with --compact-vignettes?",
        "Set vignette's cache to FALSE?"
    )
}
