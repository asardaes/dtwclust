#' Time series clustering along with optimizations for the Dynamic Time Warping distance
#'
#' Time series clustering with a wide variety of strategies and a series of optimizations specific to
#' the Dynamic Time Warping (DTW) distance and its corresponding lower bounds (LBs). There are
#' implementations of both traditional clustering algorithms, and more recent procedures such as
#' k-Shape and TADPole clustering. The functionality can be easily extended with custom distance
#' measures and centroid definitions.
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
#' Most optimizations require equal dimensionality, which means time series should have equal length.
#' DTW itself does not require this, but it is relatively expensive to compute. Other distance definitions
#' may be used, or series could be reinterpolated to a matching length (Ratanamahatana and Keogh,
#' 2004).
#'
#' Other packages that are particularly leveraged here are the \code{proxy} package for distance
#' matrix calculations, and the \code{dtw} package for the core DTW calculations. The main clustering
#' function and entry point for this package is \code{\link{dtwclust}}. There are a couple of things
#' to bear in mind while using this package:
#'
#' Most distance calculations make use of the \code{\link[proxy]{dist}} function in the \code{proxy}
#' package, which accepts one or two arguments for data objects (\code{x} and \code{y}). Users should
#' use the two-input \strong{list} version, even if there is just one dataset (i.e.
#' \code{proxy::dist(x = dataset, y = dataset, ...)}), because otherwise it sometimes fails to
#' detect a whole time series as a single object and, instead, calculates distances between each observation
#' of each time series.
#'
#' Additionally, please note the random number generator is set to L'Ecuyer-CMRG when \code{dtwclust}
#' is attached in an attempt to preserve reproducibility. You are free to change this afterwards if
#' you wish. See \code{\link[base]{RNGkind}}.
#'
#' @docType package
#' @name dtwclust-package
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Begum N, Ulanova L, Wang J and Keogh E (2015). ``Accelerating Dynamic Time Warping Clustering with
#' a Novel Admissible Pruning Strategy.'' In \emph{Conference on Knowledge Discovery and Data Mining},
#' series KDD '15. ISBN 978-1-4503-3664-2/15/08, \url{http://dx.doi.org/10.1145/2783258.2783286}.
#'
#' Giorgino T (2009). ``Computing and Visualizing Dynamic Time Warping Alignments in \code{R}: The
#' 'dtw' Package.'' \emph{Journal of Statistical Software}, \strong{31}(7), pp. 1-24.
#' \url{http://www.jstatsoft.org/v31/i07/}.
#'
#' Ratanamahatana A and Keogh E (2004). ``Everything you know about dynamic time warping is wrong.''
#' In \emph{3rd Workshop on Mining Temporal and Sequential Data, in conjunction with 10th ACM SIGKDD
#' Int. Conf. Knowledge Discovery and Data Mining (KDD-2004), Seattle, WA}.
#'
#' Keogh E and Ratanamahatana CA (2005). ``Exact indexing of dynamic time warping.'' \emph{Knowledge
#' and information systems}, \strong{7}(3), pp. 358-386.
#'
#' Lemire D (2009). ``Faster retrieval with a two-pass dynamic-time-warping lower bound .''
#' \emph{Pattern Recognition}, \strong{42}(9), pp. 2169 - 2180. ISSN 0031-3203,
#' \url{http://dx.doi.org/10.1016/j.patcog.2008.11.030},
#' \url{http://www.sciencedirect.com/science/article/pii/S0031320308004925}.
#'
#' Liao TW (2005). ``Clustering of time series data - a survey.'' \emph{Pattern recognition},
#' \strong{38}(11), pp. 1857-1874.
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In \emph{Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data},
#' series SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' Petitjean F, Ketterlin A and Gancarski P (2011). ``A global averaging method for dynamic time
#' warping, with applications to clustering.'' \emph{Pattern Recognition}, \strong{44}(3), pp. 678 -
#' 693. ISSN 0031-3203, \url{http://dx.doi.org/10.1016/j.patcog.2010.09.013},
#' \url{http://www.sciencedirect.com/science/article/pii/S003132031000453X}.
#'
#' Sakoe H and Chiba S (1978). ``Dynamic programming algorithm optimization for spoken word
#' recognition.'' \emph{Acoustics, Speech and Signal Processing, IEEE Transactions on},
#' \strong{26}(1), pp. 43-49. ISSN 0096-3518, \url{http://doi.org/10.1109/TASSP.1978.1163055}.
#'
#' Wang X, Mueen A, Ding H, Trajcevski G, Scheuermann P and Keogh E (2013). ``Experimental comparison
#' of representation methods and distance measures for time series data.'' \emph{Data Mining and
#' Knowledge Discovery}, \strong{26}(2), pp. 275-309. ISSN 1384-5810,
#' \url{http://doi.org/10.1007/s10618-012-0250-5}, \url{http://dx.doi.org/10.1007/s10618-012-0250-5}.
#'
#' Bedzek, J.C. (1981). Pattern recognition with fuzzy objective function algorithms.
#'
#' D'Urso, P., & Maharaj, E. A. (2009). Autocorrelation-based fuzzy clustering of time series.
#' Fuzzy Sets and Systems, 160(24), 3565-3589.
#'
#' @seealso
#'
#' \code{\link{dtwclust}}, \code{\link[proxy]{dist}}, \code{\link[dtw]{dtw}}
#'
#' @include utils.R
#'
#' @useDynLib dtwclust
#'
#' @import methods
#' @import proxy
#' @import foreach
#' @import ggplot2
#'
#' @importFrom dtw dtw
#' @importFrom reshape2 melt
#' @importFrom parallel splitIndices
#' @importFrom caTools runmin
#' @importFrom caTools runmax
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
#' @importFrom flexclust randIndex
#' @importFrom flexclust clusterSim
#' @importFrom graphics plot
#' @importFrom rngtools RNGseq
#' @importFrom rngtools setRNG
#' @importFrom utils head
#'
NULL

.onAttach <- function(lib, pkg) {

     ## Register DTW2
     if (!consistency_check("DTW2", "dist", silent = TRUE))
          proxy::pr_DB$set_entry(FUN = dtw2, names=c("DTW2", "dtw2"),
                                 loop = TRUE, type = "metric", distance = TRUE,
                                 description = "DTW with L2 norm",
                                 PACKAGE = "dtwclust")

     ## Register LB_Keogh with the 'proxy' package for distance matrix calculation
     if (!consistency_check("LB_Keogh", "dist", silent = TRUE))
          proxy::pr_DB$set_entry(FUN = lb_keogh_loop, names=c("LBK", "LB_Keogh", "lbk"),
                                 loop = FALSE, type = "metric", distance = TRUE,
                                 description = "Keogh's DTW lower bound for the Sakoe-Chiba band",
                                 PACKAGE = "dtwclust") #, PREFUN = proxy_prefun)


     ## Register LB_Improved with the 'proxy' package for distance matrix calculation
     if (!consistency_check("LB_Improved", "dist", silent = TRUE))
          proxy::pr_DB$set_entry(FUN = lb_improved_loop, names=c("LBI", "LB_Improved", "lbi"),
                                 loop = FALSE, type = "metric", distance = TRUE,
                                 description = "Lemire's improved DTW lower bound for the Sakoe-Chiba band",
                                 PACKAGE = "dtwclust") #, PREFUN = proxy_prefun)

     ## Register SBD
     if (!consistency_check("SBD", "dist", silent = TRUE))
          proxy::pr_DB$set_entry(FUN = SBD.proxy, names=c("SBD", "sbd"),
                                 loop = FALSE, type = "metric", distance = TRUE,
                                 description = "Paparrizos and Gravanos' shape-based distance for time series",
                                 PACKAGE = "dtwclust") #, PREFUN = proxy_prefun)

     ## Register DTW_LB
     if (!consistency_check("DTW_LB", "dist", silent = TRUE))
          proxy::pr_DB$set_entry(FUN = dtw_lb, names=c("DTW_LB", "dtw_lb"),
                                 loop = FALSE, type = "metric", distance = TRUE,
                                 description = "DTW distance aided with Lemire's lower bound",
                                 PACKAGE = "dtwclust") #, PREFUN = proxy_prefun)

     RNGkind("L'Ecuyer")

     packageStartupMessage("Setting random number generator to L'Ecuyer-CMRG")
}

.onUnload <- function(libpath) {
     library.dynam.unload("dtwclust", libpath)
}
