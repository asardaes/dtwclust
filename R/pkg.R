#' Time series clustering under Dynamic Time Warping (DTW) distance.
#'
#' Perform time series clustering using different techniques related to the DTW distance and its
#' corresponding lower bounds (LB). Additionally, an implementation of k-Shape clustering is available.
#'
#' This package tries to consolidate the different procedures available to perform clustering of time series
#' under DTW. Right now only univariate time series are supported. Similarly, time series should have equal
#' lengths for partitional methods.
#'
#' Please see the documentation for \code{\link{dtwclust}}, which serves as the main entry point.
#'
#' Other packages that are particularly leveraged here are the \code{flexclust} package for partitional
#' clustering, the \code{proxy} package for distance matrix calculations, and the \code{dtw} package for the
#' core DTW calculations.
#'
#' Four distances are registered via \code{\link[proxy]{pr_DB}}: \code{"LB_Keogh", "LB_Improved", "SBD"} and
#' \code{"DTW2"}. See \code{\link{lb_keogh}}, \code{\link{lb_improved}} and \code{\link{SBD}} for more
#' details on the first 3. The last one is done with \code{\link[dtw]{dtw}} using \code{L2} norm, but it
#' differs from the result you would obtain if you specify \code{L2} as \code{dist.method}: with \code{DTW2},
#' pointwise distances (the local cost matrix) are calculated with \code{L1} norm, \emph{each} element of the
#' matrix is squared and the result is fed into \code{\link[dtw]{dtw}}, which finds the optimum warping path.
#' The square root of the resulting distance is \emph{then} computed.
#'
#' Please note that the \code{\link[proxy]{dist}} function in the \code{proxy} package accepts one or two
#' arguments for data objects. Users should usually use the two-input version, even if there is just one
#' dataset (i.e. \code{proxy::dist(x=data, y=data, ...)}), because the one-input version sometimes fails to
#' detect a whole time series as a single object and, instead, calculates distances between each observation
#' of each time series.
#'
#' @docType package
#' @name dtwclust-package
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Begum N, Ulanova L, Wang J and Keogh E (2015). ``Accelerating Dynamic Time Warping Clustering with a Novel Admissible Pruning
#' Strategy.'' In \emph{Conference on Knowledge Discovery and Data Mining}, series KDD '15. ISBN 978-1-4503-3664-2/15/08, \url{
#' http://www.cs.ucr.edu/~eamonn/Speeded\%20Clustering\%20Paper\%20Camera\%20Ready.pdf}.
#'
#' Giorgino T (2009). ``Computing and Visualizing Dynamic Time Warping Alignments in \code{R}: The 'dtw' Package.'' \emph{Journal
#' of Statistical Software}, \strong{31}(7), pp. 1-24. \url{http://www.jstatsoft.org/v31/i07/}.
#'
#' Ratanamahatana A and Keogh E (2004). ``Everything you know about dynamic time warping is wrong.'' In \emph{3rd Workshop on Mining
#' Temporal and Sequential Data, in conjunction with 10th ACM SIGKDD Int. Conf. Knowledge Discovery and Data Mining (KDD-2004),
#' Seattle, WA}.
#'
#' Keogh E and Ratanamahatana CA (2005). ``Exact indexing of dynamic time warping.'' \emph{Knowledge and information systems}, \strong{7}(3),
#' pp. 358-386.
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .'' \emph{Pattern Recognition}, \strong{42}(9), pp.
#' 2169 - 2180. ISSN 0031-3203, \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030}, \url{
#' http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
#'
#' Liao TW (2005). ``Clustering of time series data - a survey.'' \emph{Pattern recognition}, \strong{38}(11), pp. 1857-1874.
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.'' In \emph{Proceedings of the 2015
#' ACM SIGMOD International Conference on Management of Data}, series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{
#' http://doi.org/10.1145/2723372.2737793}.
#'
#' Petitjean F, Ketterlin A and Gancarski P (2011). ``A global averaging method for dynamic time warping, with applications to
#' clustering.'' \emph{Pattern Recognition}, \strong{44}(3), pp. 678 - 693. ISSN 0031-3203, \url{
#' http://dx.doi.org/10.1016/j.patcog.2010.09.013}, \url{
#' http://www.sciencedirect.com/science/article/pii/S003132031000453X}.
#'
#' Sakoe H and Chiba S (1978). ``Dynamic programming algorithm optimization for spoken word recognition.'' \emph{Acoustics, Speech and
#' Signal Processing, IEEE Transactions on}, \strong{26}(1), pp. 43-49. ISSN 0096-3518, \url{
#' http://doi.org/10.1109/TASSP.1978.1163055}.
#'
#' Wang X, Mueen A, Ding H, Trajcevski G, Scheuermann P and Keogh E (2013). ``Experimental comparison of representation methods
#' and distance measures for time series data.'' \emph{Data Mining and Knowledge Discovery}, \strong{26}(2), pp. 275-309. ISSN 1384-5810,
#' \url{http://doi.org/10.1007/s10618-012-0250-5}, \url{http://dx.doi.org/10.1007/s10618-012-0250-5}.
#'
#' @seealso \code{\link{dtwclust}}, \code{\link[flexclust]{kcca}}, \code{\link[proxy]{dist}},
#' \code{\link[dtw]{dtw}}
NULL

.onAttach <- function(lib, pkg) {

     ## Register DTW2

     if (proxy::pr_DB$entry_exists("DTW2"))
          proxy::pr_DB$delete_entry("DTW2")

     proxy::pr_DB$set_entry(FUN = dtw2, names=c("DTW2"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "DTW with L2 as pointwise norm")

     ## Register LB_Keogh with the 'proxy' package for distance matrix calculation

     if (proxy::pr_DB$entry_exists("LB_Keogh"))
          proxy::pr_DB$delete_entry("LB_keogh")

     proxy::pr_DB$set_entry(FUN = lb_keogh_loop, names=c("LBK", "LB_Keogh"),
                            loop = FALSE, type = "metric", distance = TRUE,
                            description = "Keogh's DTW lower bound but using L1 norm")


     ## Register LB_Improved with the 'proxy' package for distance matrix calculation

     if (proxy::pr_DB$entry_exists("LB_Improved"))
          proxy::pr_DB$delete_entry("LB_Improved")

     proxy::pr_DB$set_entry(FUN = lb_improved_loop, names=c("LBI", "LB_Improved"),
                            loop = FALSE, type = "metric", distance = TRUE,
                            description = "Lemire's improved DTW lower bound using L1 norm")

     ## Register SBD

     if (proxy::pr_DB$entry_exists("SBD"))
          proxy::pr_DB$delete_entry("SBD")

     proxy::pr_DB$set_entry(FUN = SBD.proxy, names=c("SBD", "shape"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "Paparrizos' shape-based distance for time series")
}
