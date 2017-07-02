#' Time series clustering along with optimizations for the Dynamic Time Warping distance
#'
#' Time series clustering with a wide variety of strategies and a series of optimizations specific
#' to the Dynamic Time Warping (DTW) distance and its corresponding lower bounds (LBs).
#'
#' @docType package
#' @name dtwclust-package
#' @include utils.R
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
#' wrapper for multiple tests in [compare_clusterings()].
#'
#' Please note the random number generator is set to L'Ecuyer-CMRG when \pkg{dtwclust} is attached
#' in an attempt to preserve reproducibility. You are free to change this afterwards if you wish.
#' See [base::RNGkind()].
#'
#' For more information, please read the included package vignette, which can be accessed by typing
#' `vignette("dtwclust")`.
#'
#' @note
#'
#' This software package was developed independently of any organization or institution that is or
#' has been associated with the author.
#'
#' There are a couple of reasons why this package will **not** work if it is not attached (via
#' [base::library()]). For the specifics, see below.
#'
#' The \pkg{methods} [package][methods::methods-package] must be attached in order for some internal
#' inheritance to work properly. This is usually done automatically by `R`, with [utils::Rscript()]
#' being an exception. This is why this package is a dependency of \pkg{dtwclust}. Additionally, the
#' included distance functions are registered with \pkg{proxy} during attachment.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references.
#'
#' @seealso
#'
#' [tsclust()], [compare_clusterings()], [proxy::dist()], [dtw::dtw()]
#'
#' @useDynLib dtwclust, .registration = TRUE
#'
#' @import foreach
#'
#' @importFrom bigmemory attach.big.matrix
#' @importFrom bigmemory big.matrix
#' @importFrom bigmemory describe
#' @importFrom bigmemory is.big.matrix
#'
#' @importFrom clue as.cl_class_ids
#' @importFrom clue as.cl_membership
#' @importFrom clue cl_class_ids
#' @importFrom clue cl_membership
#' @importFrom clue is.cl_dendrogram
#' @importFrom clue is.cl_hard_partition
#' @importFrom clue is.cl_hierarchy
#' @importFrom clue is.cl_partition
#' @importFrom clue n_of_classes
#' @importFrom clue n_of_objects
#'
#' @importFrom dtw dtw
#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#'
#' @importFrom flexclust clusterSim
#' @importFrom flexclust comPart
#' @importFrom flexclust randIndex
#'
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#'
#' @importFrom graphics plot
#'
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom Matrix forceSymmetric
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#'
#' @importFrom methods S3Part
#' @importFrom methods as
#' @importFrom methods callNextMethod
#' @importFrom methods initialize
#' @importFrom methods new
#' @importFrom methods setAs
#' @importFrom methods setClass
#' @importFrom methods setClassUnion
#' @importFrom methods setGeneric
#' @importFrom methods setValidity
#' @importFrom methods show
#' @importFrom methods signature
#' @importFrom methods slot
#' @importFrom methods slot<-
#' @importFrom methods slotNames
#' @importFrom methods validObject
#'
#' @importFrom parallel splitIndices
#'
#' @importFrom plyr rbind.fill
#'
#' @importFrom proxy dist
#' @importFrom proxy pr_DB
#'
#' @importFrom RSpectra eigs_sym
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
#' @importFrom stats as.dist
#' @importFrom stats as.hclust
#' @importFrom stats convolve
#' @importFrom stats cutree
#' @importFrom stats fft
#' @importFrom stats hclust
#' @importFrom stats median
#' @importFrom stats nextn
#' @importFrom stats predict
#' @importFrom stats runif
#' @importFrom stats update
#'
#' @importFrom utils head
#' @importFrom utils packageVersion
#'
NULL ## remember to check methods imports after removing dtwclust()

#' Deprecated functionality in \pkg{dtwclust}
#'
#' @name dtwclust-deprecated
#' @aliases dtwclust-deprecated
#'
#' @description
#'
#' The following functions and/or classes are deprecated and will be eventually removed:
#'
#' - [dtwclust()] - see [tsclust()]
#' - [dtwclustFamily-class] - see [tsclustFamily-class]
#' - [dtwclustControl-class] - see [tsclust-controls]
#' - [dtwclust-class] - see [TSClusters-class]
#' - [dtwclust-methods] - see [tsclusters-methods]
#' - [randIndex()] - see [cvi()]
#' - [clusterSim()]
#' - [create_dtwclust()] - see [tsclusters-methods]
#' - [compute_envelop()] - see [compute_envelope()]
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

    ## Register LB_Keogh
    if (!check_consistency("LB_Keogh", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_keogh_proxy, names=c("LBK", "LB_Keogh", "lbk"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Keogh's DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)


    ## Register LB_Improved
    if (!check_consistency("LB_Improved", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_improved_proxy, names=c("LBI", "LB_Improved", "lbi"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Lemire's improved DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    ## Register SBD
    if (!check_consistency("SBD", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = SBD_proxy, names=c("SBD", "sbd"),
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
                               description = "Fast (triangular) global alignment kernel distance",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 1 - d })

    ## Register uGAK
    if (!check_consistency("uGAK", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = GAK_simil, names=c("uGAK", "ugak"),
                               loop = FALSE, type = "metric", distance = FALSE,
                               description = "Fast (triangular) global alignment kernel similarity",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    RNGkind("L'Ecuyer")

    ## avoids default message if no backend exists
    if (is.null(foreach::getDoParName())) foreach::registerDoSEQ()

    packageStartupMessage("\ndtwclust:\n",
                          "Setting random number generator to L'Ecuyer-CMRG (see RNGkind()).\n",
                          'To read the included vignette type: vignette("dtwclust").\n',
                          'Please see news(package = "dtwclust") for important information.')

    if (grepl("\\.9000$", utils::packageVersion("dtwclust")))
        packageStartupMessage("\nThis is a developer version of 'dtwclust'.\n",
                              "Using devtools::test() is currently broken, see tests/testthat.R")
}

.onUnload <- function(libpath) { library.dynam.unload("dtwclust", libpath) }

release_questions <- function() {
    c(
        "Changed .Rbuildignore to exclude test rds files?",
        "Built the binary with --compact-vignettes?",
        "Set vignette's cache to FALSE?"
    )
}
