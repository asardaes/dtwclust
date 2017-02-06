#' Time series clustering along with optimizations for the Dynamic Time Warping distance
#'
#' Time series clustering with a wide variety of strategies and a series of optimizations specific
#' to the Dynamic Time Warping (DTW) distance and its corresponding lower bounds (LBs). There are
#' implementations of both traditional clustering algorithms, and more recent procedures such as
#' k-Shape and TADPole clustering. Functionality can be easily extended with custom distance
#' measures and centroid definitions.
#'
#' @docType package
#' @name dtwclust-package
#' @include utils.R
#'
#' @details
#'
#' Many of the algorithms implemented in this package are specifically tailored to time series and
#' DTW, hence its name. However, the main clustering function is flexible so that one can test many
#' different clustering approaches, using either the time series directly, or by applying suitable
#' transformations and then clustering in the resulting space.
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
#' Other packages that are particularly leveraged here are the \pkg{proxy} package for distance
#' matrix calculations and the \pkg{dtw} package for some of the core DTW calculations.
#'
#' The main clustering function and entry point for this package is \code{\link{tsclust}}.
#'
#' Please note the random number generator is set to L'Ecuyer-CMRG when \pkg{dtwclust} is attached
#' in an attempt to preserve reproducibility. You are free to change this afterwards if you wish.
#' See \code{\link[base]{RNGkind}}.
#'
#' For more information, please read the included package vignette, which can be accessed by typing
#' \code{vignette("dtwclust")}.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references.
#'
#' @seealso
#'
#' \code{\link{tsclust}}, \code{\link[proxy]{dist}}, \code{\link[dtw]{dtw}}
#'
#' @useDynLib dtwclust
#'
#' @import clue
#' @import foreach
#' @import ggplot2
#' @import methods
#'
#' @importFrom dtw dtw
#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#'
#' @importFrom caTools runmin
#' @importFrom caTools runmax
#'
#' @importFrom flexclust randIndex
#' @importFrom flexclust clusterSim
#' @importFrom flexclust comPart
#'
#' @importFrom graphics plot
#'
#' @importFrom parallel splitIndices
#'
#' @importFrom proxy dist
#' @importFrom proxy pr_DB
#'
#' @importFrom Rcpp evalCpp
#'
#' @importFrom reshape2 melt
#'
#' @importFrom rngtools RNGseq
#' @importFrom rngtools setRNG
#'
#' @importFrom stats aggregate
#' @importFrom stats approx
#' @importFrom stats convolve
#' @importFrom stats cutree
#' @importFrom stats fft
#' @importFrom stats hclust
#' @importFrom stats median
#' @importFrom stats nextn
#' @importFrom stats update
#' @importFrom stats predict
#' @importFrom stats runif
#' @importFrom stats as.dist
#' @importFrom stats as.hclust
#'
#' @importFrom utils packageVersion
#' @importFrom utils head
#'
NULL

.onAttach <- function(lib, pkg) {
    ## proxy_prefun is in utils.R

    ## Register DTW2
    if (!check_consistency("DTW2", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw2.proxy, names=c("DTW2", "dtw2"),
                               loop = TRUE, type = "metric", distance = TRUE,
                               description = "DTW with L2 norm",
                               PACKAGE = "dtwclust")

    ## Register DTW_BASIC
    if (!check_consistency("DTW_BASIC", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw_basic_proxy, names=c("DTW_BASIC", "dtw_basic"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Basic and maybe faster DTW distance",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    ## Register LB_Keogh with the 'proxy' package for distance matrix calculation
    if (!check_consistency("LB_Keogh", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_keogh_proxy, names=c("LBK", "LB_Keogh", "lbk"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Keogh's DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)


    ## Register LB_Improved with the 'proxy' package for distance matrix calculation
    if (!check_consistency("LB_Improved", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_improved_proxy, names=c("LBI", "LB_Improved", "lbi"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Lemire's improved DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    ## Register SBD
    if (!check_consistency("SBD", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = SBD.proxy, names=c("SBD", "sbd"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Paparrizos and Gravanos' shape-based distance for time series",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 2 - d })

    ## Register DTW_LB
    if (!check_consistency("DTW_LB", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw_lb, names=c("DTW_LB", "dtw_lb"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "DTW distance aided with Lemire's lower bound",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    ## Register GAK
    if (!check_consistency("GAK", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = GAK_proxy, names=c("GAK", "gak"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Fast (triangular) global alignment kernel",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 1 - d })

    RNGkind("L'Ecuyer")

    ## avoids default message if no backend exists
    if (is.null(foreach::getDoParName()))
        foreach::registerDoSEQ()

    packageStartupMessage("\ndtwclust: Setting random number generator to L'Ecuyer-CMRG (see RNGkind()).\n",
                          'To read the included vignette, type: vignette("dtwclust").\n',
                          'Please see news(package = "dtwclust") for important information!\n')

    if (grepl("\\.9000$", utils::packageVersion("dtwclust")))
        packageStartupMessage("This is a developer version of 'dtwclust'.")
}

.onUnload <- function(libpath) {
    library.dynam.unload("dtwclust", libpath)
}

release_questions <- function() {
    c(
        "Changed .Rbuildignore to exclude test rds files?",
        "Built the binary with --compact-vignettes?"
    )
}
