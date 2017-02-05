#' Cluster validity indices
#'
#' Compute different cluster validity indices (CVIs) of a given cluster partition, using the
#' clustering distance measure and centroid function if applicable.
#'
#' @export
#' @exportMethod cvi
#'
#' @param a An object returned by the \code{\link{dtwclust}} or \code{\link{tsclust}} function, or a
#'   vector that can be coerced to integers which indicate the cluster memeberships.
#' @param b If needed, a vector that can be coerced to integers which indicate the cluster
#'   memeberships. The ground truth (if known) should be provided here.
#' @param type Character vector indicating which indices are to be computed. See supported values
#'   below.
#' @param ... Arguments to pass to and from other methods.
#' @param log.base Base of the logarithm to be used in the calculation of VI.
#'
#' @details
#'
#' Clustering is commonly considered to be an unsupervised procedure, so evaluating its performance
#' can be rather subjective. However, a great amount of effort has been invested in trying to
#' standardize cluster evaluation metrics by using cluster validity indices (CVIs).
#'
#' CVIs can be classified as internal, external or relative depending on how they are computed.
#' Focusing on the first two, the crucial difference is that internal CVIs only consider the
#' partitioned data and try to define a measure of cluster purity, whereas external CVIs compare the
#' obtained partition to the correct one. Thus, external CVIs can only be used if the ground truth
#' is known. Each index defines their range of values and whether they are to be minimized or
#' maximized. In many cases, these CVIs can be used to evaluate the result of a clustering algorithm
#' regardless of how the clustering works internally, or how the partition came to be.
#'
#' Knowing which CVI will work best cannot be determined a priori, so they should be tested for each
#' specific application. Usually, many CVIs are utilized and compared to each other, maybe using a
#' majority vote to decide on a final result. Furthermore, it should be noted that many CVIs perform
#' additional distance calculations when being computed, which can be very considerable if using
#' DTW.
#'
#' Note that, even though a fuzzy partition can be changed into a crisp one, making it compatible
#' with many of the existing CVIs, there are also fuzzy CVIs tailored specifically to fuzzy
#' clustering, and these may be more suitable in those situations, but have not been implemented
#' here yet.
#'
#' @return The chosen CVIs
#'
#' @section External CVIs:
#'
#'   The first 4 CVIs are calculated via \code{\link[flexclust]{comPart}}, so please refer to that
#'   function.
#'
#'   \itemize{
#'     \item \code{"RI"}: Rand Index (to be maximized).
#'     \item \code{"ARI"}: Adjusted Rand Index (to be maximized).
#'     \item \code{"J"}: Jaccard Index (to be maximized).
#'     \item \code{"FM"}: Fowlkes-Mallows (to be maximized).
#'     \item \code{"VI"}: Variation of Information (Meila (2003); to be minimized).
#'   }
#'
#' @section Internal CVIs:
#'
#'   The indices marked with an exclamation mark (!) calculate (or re-use if already available) the
#'   whole distance matrix between the series in the data. If you were trying to avoid this in the
#'   first place, then these CVIs might not be suitable for your application.
#'
#'   The indices marked with a question mark (?) depend on the extracted centroids, so bear that in
#'   mind if a hierarchical procedure was used and/or the centroid function has associated
#'   randomness (such as \code{\link{shape_extraction}} with series of different length).
#'
#'   The indices marked with a tilde (~) require the calculation of a global centroid. Since
#'   \code{\link{DBA}} and \code{\link{shape_extraction}} (for series of different length) have some
#'   randomness associated, these indices might not be appropriate for those centroids.
#'
#'   \itemize{
#'     \item \code{"Sil"} (!): Silhouette index (Arbelaitz et al. (2013); to be maximized).
#'     \item \code{"D"} (!): Dunn index (Arbelaitz et al. (2013); to be maximized).
#'     \item \code{"COP"} (!): COP index (Arbelaitz et al. (2013); to be minimized).
#'     \item \code{"DB"} (?): Davies-Bouldin index (Arbelaitz et al. (2013); to be minimized).
#'     \item \code{"DBstar"} (?): Modified Davies-Bouldin index (DB*) (Kim and Ramakrishna (2005);
#'     to be minimized).
#'     \item \code{"CH"} (~): Calinski-Harabasz index (Arbelaitz et al. (2013); to be maximized).
#'     \item \code{"SF"} (~): Score Function (Saitta et al. (2007); to be maximized).
#'   }
#'
#' @section Additionally:
#'
#'   \itemize{
#'     \item \code{"valid"}: Returns all valid indices depending on the type of \code{a} and whether
#'       \code{b} was provided or not.
#'     \item \code{"internal"}: Returns all internal CVIs. Only supported for
#'       \code{\link{dtwclust-class}} objects.
#'     \item \code{"external"}: Returns all external CVIs. Requires \code{b} to be provided.
#'   }
#'
#' @note
#'
#' In the original definition of many internal CVIs, the Euclidean distance and a mean centroid was
#' used. The implementations here change this, making use of whatever distance/centroid was chosen
#' during clustering.
#'
#' Some internal indices require the original data for calculations, so the control flag
#' \code{save.data} must be set to \code{TRUE} when running the clustering algorithm.
#'
#' @references
#'
#' Arbelaitz, O., Gurrutxaga, I., Muguerza, J., Perez, J. M., & Perona, I. (2013). An extensive
#' comparative study of cluster validity indices. Pattern Recognition, 46(1), 243-256.
#'
#' Kim, M., & Ramakrishna, R. S. (2005). New indices for cluster validity assessment. Pattern
#' Recognition Letters, 26(15), 2353-2363.
#'
#' Meila, M. (2003). Comparing clusterings by the variation of information. In Learning theory and
#' kernel machines (pp. 173-187). Springer Berlin Heidelberg.
#'
#' Saitta, S., Raphael, B., & Smith, I. F. (2007). A bounded index for cluster validity. In
#' International Workshop on Machine Learning and Data Mining in Pattern Recognition (pp. 174-187).
#' Springer Berlin Heidelberg.
#'
setGeneric("cvi", def = function(a, b = NULL, type = "valid", ..., log.base = 10) {
    ## Only external CVIs for default, dtwclust-methods.R has the internal ones
    if (is.null(b))
        stop("A second set of cluster membership indices is required in 'b' for this/these CVI(s).")

    a <- as.integer(a)
    b <- as.integer(b)

    if (length(a) != length(b))
        stop("External CVIs: the length of 'a' and 'b' must match.")

    type <- match.arg(type, several.ok = TRUE,
                      c("RI", "ARI", "J", "FM", "VI",
                        "valid", "external"))

    if (any(type %in% c("valid", "external")))
        type <- c("RI", "ARI", "J", "FM", "VI")

    which_flexclust <- type %in% c("RI", "ARI", "J", "FM")

    if (any(which_flexclust))
        CVIs <- flexclust::comPart(x = a, y = b, type = type[which_flexclust])
    else
        CVIs <- numeric()

    if (any(type == "VI")) {
        ## Variation of information
        ## taken from https://github.com/cran/mcclust/blob/master/R/vi.dist.R

        ## entropy
        ent <- function(cl) {
            n <- length(cl)
            p <- table(cl) / n
            -sum(p * log(p, base = log.base))
        }

        ## mutual information
        mi <- function(cl1, cl2) {
            p12 <- table(cl1, cl2) / length(cl1)
            p1p2 <- outer(table(cl1) / length(cl1), table(cl2) / length(cl2))
            sum(p12[p12 > 0] * log(p12[p12 > 0] / p1p2[p12 > 0], base = log.base))
        }

        VI <- ent(a) + ent(b) - 2 * mi(a, b)
        CVIs <- c(CVIs, VI = VI)
    }

    CVIs
})
