#' Time series clustering with Dynamic Time Warping (DTW) distance.
#'
#' Something
#'
#' Something more
#'
#' @docType package
#' @name dtwclust-package
NULL

.onAttach <- function(lib, pkg) {

     ## Register DTW2

     if (proxy::pr_DB$entry_exists("DTW2"))
          proxy::pr_DB$delete_entry("DTW2")

     proxy::pr_DB$set_entry(FUN = dtw2, names=c("DTW2"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "DTW with L2 as pointwise norm")

     ## Register LB_Keogh with the 'proxy' package for distance matrix calculation

     if (proxy::pr_DB$entry_exists("LB_Keogh1"))
          proxy::pr_DB$delete_entry("LB_keogh1")

     proxy::pr_DB$set_entry(FUN = lb_keogh_loop, names=c("LB_Keogh1", "LBK", "LB_Keogh"),
                     loop = FALSE, type = "metric", distance = TRUE,
                     description = "Keogh's DTW lower bound but using L1 norm")


     ## Register LB_Improved with the 'proxy' package for distance matrix calculation

     if (proxy::pr_DB$entry_exists("LB_Improved1"))
          proxy::pr_DB$delete_entry("LB_Improved1")

     proxy::pr_DB$set_entry(FUN = lb_improved_loop, names=c("LB_Improved1", "LBI", "LB_Improved"),
                     loop = FALSE, type = "metric", distance = TRUE,
                     description = "Lemire's improved DTW lower bound using L1 norm")

     ## Register SBD

     if (proxy::pr_DB$entry_exists("SBD"))
          proxy::pr_DB$delete_entry("SBD")

     proxy::pr_DB$set_entry(FUN = SBD.proxy, names=c("SBD", "shape"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "Paparrizos' shape-based distance for time series")
}
