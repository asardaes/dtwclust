#' Cluster validity indices
#'
#' Compute different cluster validity indices (CVIs) of a given cluster partition, using the
#' clustering distance measure and centroid function if applicable.
#'
#' @export
#' @exportMethod cvi
#' @importFrom flexclust comPart
#' @importFrom methods setGeneric
#'
#' @param a An object returned by [tsclust()], for crisp partitions a vector that can be coerced to
#'   integers which indicate the cluster memeberships, or the membership matrix for soft clustering.
#' @param b If needed, a vector that can be coerced to integers which indicate the cluster
#'   memeberships. The ground truth (if known) should be provided here.
#' @param type Character vector indicating which indices are to be computed. See supported values
#'   below.
#' @param ... Arguments to pass to and from other methods.
#' @param log.base Base of the logarithm to be used in the calculation of VI (see details).
#'
#' @details
#'
#' Clustering is commonly considered to be an unsupervised procedure, so evaluating its performance
#' can be rather subjective. However, a great amount of effort has been invested in trying to
#' standardize cluster evaluation metrics by using cluster validity indices (CVIs).
#'
#' In general, CVIs can be either tailored to crisp or fuzzy partitions. CVIs can be classified as
#' internal, external or relative depending on how they are computed. Focusing on the first two, the
#' crucial difference is that internal CVIs only consider the partitioned data and try to define a
#' measure of cluster purity, whereas external CVIs compare the obtained partition to the correct
#' one. Thus, external CVIs can only be used if the ground truth is known.
#'
#' Note that even though a fuzzy partition can be changed into a crisp one, making it compatible
#' with many of the existing crisp CVIs, there are also fuzzy CVIs tailored specifically to fuzzy
#' clustering, and these may be more suitable in those situations. Fuzzy partitions usually have no
#' ground truth associated with them, but there are exceptions depending on the task's goal.
#'
#' Each index defines their range of values and whether they are to be minimized or maximized. In
#' many cases, these CVIs can be used to evaluate the result of a clustering algorithm regardless of
#' how the clustering works internally, or how the partition came to be.
#'
#' Knowing which CVI will work best cannot be determined a priori, so they should be tested for each
#' specific application. Usually, many CVIs are utilized and compared to each other, maybe using a
#' majority vote to decide on a final result. Furthermore, it should be noted that many CVIs perform
#' additional distance calculations when being computed, which can be very considerable if using DTW
#' or GAK.
#'
#' @return The chosen CVIs
#'
#' @section External CVIs:
#'
#'   - Crisp partitions (the first 4 are calculated via [flexclust::comPart()])
#'     + `"RI"`: Rand Index (to be maximized).
#'     + `"ARI"`: Adjusted Rand Index (to be maximized).
#'     + `"J"`: Jaccard Index (to be maximized).
#'     + `"FM"`: Fowlkes-Mallows (to be maximized).
#'     + `"VI"`: Variation of Information (Meila (2003); to be minimized).
#'   - Fuzzy partitions (based on Lei et al. (2017))
#'     + `"RI"`: Soft Rand Index (to be maximized).
#'     + `"ARI"`: Soft Adjusted Rand Index (to be maximized).
#'     + `"VI"`: Soft Variation of Information (to be minimized).
#'     + `"NMIM"`: Soft Normalized Mutual Information based on Max entropy (to be maximized).
#'
#' @section Internal CVIs:
#'
#'   The indices marked with an exclamation mark (!) calculate (or re-use if already available) the
#'   whole distance matrix between the series in the data. If you were trying to avoid this in the
#'   first place, then these CVIs might not be suitable for your application.
#'
#'   The indices marked with a question mark (?) depend on the extracted centroids, so bear that in
#'   mind if a hierarchical procedure was used and/or the centroid function has associated
#'   randomness (such as [shape_extraction()] with series of different length).
#'
#'   The indices marked with a tilde (~) require the calculation of a global centroid. Since [DBA()]
#'   and [shape_extraction()] (for series of different length) have some randomness associated,
#'   these indices might not be appropriate for those centroids.
#'
#'   - Crisp partitions
#'     + `"Sil"` (!): Silhouette index (Rousseeuw (1987); to be maximized).
#'     + `"D"` (!): Dunn index (Arbelaitz et al. (2013); to be maximized).
#'     + `"COP"` (!): COP index (Arbelaitz et al. (2013); to be minimized).
#'     + `"DB"` (?): Davies-Bouldin index (Arbelaitz et al. (2013); to be minimized).
#'     + `"DBstar"` (?): Modified Davies-Bouldin index (DB*) (Kim and Ramakrishna (2005); to be
#'       minimized).
#'     + `"CH"` (~): Calinski-Harabasz index (Arbelaitz et al. (2013); to be maximized).
#'     + `"SF"` (~): Score Function (Saitta et al. (2007); to be maximized; see notes).
#'   - Fuzzy partitions (using the nomenclature used in Wang and Zhang (2007))
#'     + `"MPC"`: to be maximized.
#'     + `"K"` (~): to be minimized.
#'     + `"T"`: to be minimized.
#'     + `"SC"` (~): to be maximized.
#'     + `"PBMF"` (~): to be maximized (see notes).
#'
#' @section Additionally:
#'
#'   - `"valid"`: Returns all valid indices depending on the type of `a` and whether `b` was
#'     provided or not.
#'   - `"internal"`: Returns all internal CVIs. Only supported for [TSClusters-class] objects.
#'   - `"external"`: Returns all external CVIs. Requires `b` to be provided.
#'
#' @note
#'
#' In the original definition of many internal and fuzzy CVIs, the Euclidean distance and a mean
#' centroid was used. The implementations here change this, making use of whatever distance/centroid
#' was chosen during clustering. However, some of the CVIs assume that the distances are symmetric,
#' since cross-distance matrices are calculated and only the upper/lower triangulars are considered.
#' A warning will be given if the matrices are not symmetric and the CVI assumes so.
#'
#' The formula for the SF index in Saitta et al. (2007) does not correspond to the one in Arbelaitz
#' et al. (2013). The one specified in the former is used here.
#'
#' The formulas for the Silhouette index are not entirely correct in Arbelaitz et al. (2013), refer
#' to Rousseeuw (1987) for the correct ones.
#'
#' The formulas for the PBMF index are not entirely unambiguous in the literature. The ones given in
#' Lin (2013) are used here.
#'
#' @references
#'
#' Arbelaitz, O., Gurrutxaga, I., Muguerza, J., Perez, J. M., & Perona, I. (2013). An extensive
#' comparative study of cluster validity indices. Pattern Recognition, 46(1), 243-256.
#'
#' Kim, M., & Ramakrishna, R. S. (2005). New indices for cluster validity assessment. Pattern
#' Recognition Letters, 26(15), 2353-2363.
#'
#' Lei, Y., Bezdek, J. C., Chan, J., Vinh, N. X., Romano, S., & Bailey, J. (2017). Extending
#' information-theoretic validity indices for fuzzy clustering. IEEE Transactions on Fuzzy Systems,
#' 25(4), 1013-1018.
#'
#' Lin, H. Y. (2013). Effective Feature Selection for Multi-class Classification Models. In
#' Proceedings of the World Congress on Engineering (Vol. 3).
#'
#' Meila, M. (2003). Comparing clusterings by the variation of information. In Learning theory and
#' kernel machines (pp. 173-187). Springer Berlin Heidelberg.
#'
#' Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of
#' cluster analysis. Journal of computational and applied mathematics, 20, 53-65.
#'
#' Saitta, S., Raphael, B., & Smith, I. F. (2007). A bounded index for cluster validity. In
#' International Workshop on Machine Learning and Data Mining in Pattern Recognition (pp. 174-187).
#' Springer Berlin Heidelberg.
#'
#' Wang, W., & Zhang, Y. (2007). On fuzzy cluster validity indices. Fuzzy sets and systems, 158(19),
#' 2095-2117.
#'
setGeneric("cvi", def = function(a, b = NULL, type = "valid", ..., log.base = 10) {
    # Only external CVIs for default, TSClusters-methods.R has the internal ones
    if (is.null(b))
        stop("A second set of cluster membership indices is required in 'b' for this/these CVI(s).")

    a <- as.integer(a)
    b <- as.integer(b)

    if (length(a) != length(b)) stop("External CVIs: the length of 'a' and 'b' must match.")

    type <- match.arg(type, several.ok = TRUE,
                      choices = c("RI", "ARI", "J", "FM", "VI", "valid", "external"))

    if (any(type %in% c("valid", "external"))) type <- c("RI", "ARI", "J", "FM", "VI")
    which_flexclust <- type %in% c("RI", "ARI", "J", "FM")

    if (any(which_flexclust))
        CVIs <- flexclust::comPart(x = a, y = b, type = type[which_flexclust])
    else
        CVIs <- numeric()

    if (any(type == "VI")) {
        # Variation of information
        # taken from https://github.com/cran/mcclust/blob/master/R/vi.dist.R

        # entropy
        ent <- function(cl) {
            p <- table(cl) / length(cl)
            -sum(p * log(p, base = log.base))
        }

        # mutual information
        mi <- function(cl1, cl2) {
            p12 <- table(cl1, cl2) / length(cl1)
            p1p2 <- outer(table(cl1) / length(cl1), table(cl2) / length(cl2))
            sum(p12[p12 > 0] * log(p12[p12 > 0] / p1p2[p12 > 0], base = log.base))
        }

        VI <- ent(a) + ent(b) - 2 * mi(a, b)
        CVIs <- c(CVIs, VI = VI)
    }
    # return
    CVIs
})

# For external fuzzy CVIs
#' @rdname cvi
#' @export
#' @exportMethod cvi
#' @importFrom methods setMethod
#' @importFrom methods signature
#'
setMethod(
    "cvi", signature = methods::signature(a = "matrix"),
    function(a, b = NULL, type = "valid", ..., log.base = 10) {
        if (is.null(b)) stop("A second set of cluster membership indices is required in 'b' for this/these CVI(s).")
        dim_b <- dim(b)
        b <- as.integer(b)
        dim(b) <- dim_b

        if (is.null(dim_b)) {
            if (nrow(a) != length(b)) stop("External CVIs: 'a'-rows and 'b'-length must match.")
            temp <- matrix(0L, nrow(a), max(b))
            temp[cbind(1L:nrow(temp), b)] <- 1L
            b <- temp
        }

        type <- match.arg(type, several.ok = TRUE,
                          choices = c("ARI", "RI", "VI", "NMIM", "valid", "external"))
        if (any(type %in% c("valid", "external"))) type <- c("ARI", "RI", "VI", "NMIM")

        num_objects <- nrow(a)
        contingency_table <- t(a) %*% b

        if (any(c("ARI", "RI") %in% type)) {
            total_sum <- sum(contingency_table)
            squared_sum <- sum(contingency_table ^ 2L)
            row_sum <- sum(base::rowSums(contingency_table) ^ 2L)
            col_sum <- sum(base::colSums(contingency_table) ^ 2L)
            pairs_in_both <- (squared_sum - total_sum) / 2
            pairs_in_neither <- (num_objects ^ 2L + squared_sum - row_sum - col_sum) / 2
            just_in_a <- (col_sum - squared_sum) / 2
            just_in_b <- (row_sum - squared_sum) / 2
        }

        if (any(c("VI", "NMIM") %in% type)) {
            joint_distribution <- contingency_table / num_objects
            joint_entropy <- -sum(joint_distribution * log(joint_distribution + .Machine$double.eps,
                                                           base = log.base))
            jdx <- base::rowSums(joint_distribution)
            jdy <- base::colSums(joint_distribution)
            entropy_x <- -sum(jdx * log(jdx + .Machine$double.eps, base = log.base))
            entropy_y <- -sum(jdy * log(jdy + .Machine$double.eps, base = log.base))
            mutual_information <- entropy_x + entropy_y - joint_entropy
        }

        CVIs <- sapply(type, function(CVI) {
            switch(EXPR = CVI,
                   # -------------------------------------------------------------------------------
                   "ARI" = {
                       sum_all <- pairs_in_both + pairs_in_neither + just_in_a + just_in_b
                       both_and_a <- pairs_in_both + just_in_a
                       both_and_b <- pairs_in_both + just_in_b
                       temp_sum <- both_and_a + both_and_b
                       temp_mul <- both_and_a * both_and_b
                       (pairs_in_both - temp_mul / sum_all) / (0.5 * temp_sum - temp_mul / sum_all)
                   },

                   # -------------------------------------------------------------------------------
                   "RI" = {
                       (pairs_in_both + pairs_in_neither) /
                           (pairs_in_both + pairs_in_neither + just_in_a + just_in_b)
                   },

                   # -------------------------------------------------------------------------------
                   "VI" = {
                       entropy_x + entropy_y - 2 * mutual_information
                   },

                   # -------------------------------------------------------------------------------
                   "NMIM" = {
                       mutual_information / max(entropy_x, entropy_y)
                   })
        })
        # return
        CVIs
    }
)
