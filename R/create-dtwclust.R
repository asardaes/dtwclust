#' Create formal `dtwclust` objects
#'
#' **Deprecated** helper function to manually create formal [dtwclust-class] objects
#'
#' @export
#' @keywords internal
#'
#' @param ... Any valid slots of [dtwclust-class].
#' @param override.family Attempt to substitute the default family with one that conforms to the
#'   provided elements? See details.
#'
#' @details
#'
#' Please see the initialize section in [tsclusters-methods] for the new objects..
#'
#' @return A [dtwclust-class] object.
#'
create_dtwclust <- function(..., override.family = TRUE) {
    tic <- proc.time()

    dots <- list(...)

    ## even if it's NULL, it'll be converted correctly
    dots$control <- as(dots$control, "dtwclustControl")

    ## some minor checks
    if (!is.null(dots$datalist)) dots$datalist <- tslist(dots$datalist)
    if (!is.null(dots$centroids)) dots$centroids <- tslist(dots$centroids)

    ## avoid infinite recursion
    if (is.null(dots$call)) {
        call <- match.call()

    } else {
        call <- dots$call
        dots$call <- NULL
    }

    .Object <- do.call(methods::new, enlist(Class = "dtwclust", dots = dots))
    .Object@call <- call

    ## some "defaults"
    if (is.null(dots$preproc)) .Object@preproc <- "none"
    if (is.null(dots$iter)) .Object@iter <- 1L
    if (is.null(dots$converged)) .Object@converged <- TRUE
    if (is.null(dots$k)) .Object@k <- length(.Object@centroids)

    ## more helpful for hierarchical/tadpole
    if (override.family) {
        if (length(.Object@type) == 0L)
            warning("Could not override family, 'type' slot is missing.")
        else if (length(.Object@distance) == 0L)
            warning("Could not override family, 'distance' slot is missing.")
        else {
            centroids <- .Object@centroids
            datalist <- .Object@datalist

            if (.Object@type == "partitional" && length(.Object@centroid))
                allcent <- all_cent(.Object@centroid,
                                    distmat = .Object@distmat,
                                    control = .Object@control,
                                    fuzzy = isTRUE(.Object@type == "fuzzy"))
            else if (.Object@type == "hierarchical" && length(formals(.Object@family@allcent)))
                allcent <- .Object@family@allcent
            else if (.Object@type == "hierarchical" && length(centroids))
                allcent <- function(dummy) {
                    datalist[which.min(apply(.Object@distmat, 1L, sum))] # for CVI's global_cent
                }
            else if (.Object@type == "tadpole" && length(centroids))
                allcent <- function(dummy) { centroids[1L] } # for CVI's global_cent
            else if (.Object@type == "fuzzy")
                allcent <- .Object@centroid
            else
                allcent <- .Object@family@allcent

            .Object@family <- new("dtwclustFamily",
                                  dist = .Object@distance,
                                  allcent = allcent,
                                  preproc = .Object@family@preproc,
                                  distmat = NULL,
                                  control = .Object@control,
                                  fuzzy = isTRUE(.Object@type == "fuzzy"))

            assign("distfun", .Object@family@dist, environment(.Object@family@allcent))

            if (.Object@type == "partitional" && .Object@centroid == "shape") {
                .Object@family@preproc <- zscore
                .Object@preproc <- "zscore"
            }
        }
    }

    if (!nrow(.Object@cldist) && length(formals(.Object@family@dist)) && length(.Object@cluster)) {
        ## no cldist available, but dist and cluster can be used to calculate it
        dm <- do.call(.Object@family@dist,
                      enlist(.Object@datalist,
                             .Object@centroids,
                             dots = .Object@dots))

        .Object@cldist <- base::as.matrix(dm[cbind(1L:length(.Object@datalist), .Object@cluster)])

        dimnames(.Object@cldist) <- NULL
    }

    if (!nrow(.Object@clusinfo) && length(.Object@cluster) && nrow(.Object@cldist)) {
        ## no clusinfo available, but cluster and cldist can be used to calculate it
        size <- as.vector(table(.Object@cluster))
        clusinfo <- data.frame(size = size, av_dist = 0)
        clusinfo[clusinfo$size > 0L, "av_dist"] <-
            as.vector(tapply(.Object@cldist[ , 1L], .Object@cluster, mean))

        .Object@clusinfo <- clusinfo
    }

    if (.Object@type == "fuzzy" && !nrow(.Object@fcluster) && length(formals(.Object@family@dist))) {
        ## no fcluster available, but dist and cluster function can be used to calculate it
        dm <- do.call(.Object@family@dist,
                      enlist(.Object@datalist,
                             .Object@centroids,
                             dots = .Object@dots))

        .Object@fcluster <- .Object@family@cluster(dm, m = .Object@control@fuzziness)
        colnames(.Object@fcluster) <- paste0("cluster_", 1:.Object@k)
    }

    ## default for when it doesn't apply
    if (.Object@type != "fuzzy") .Object@fcluster <- matrix(NA_real_)

    ## just a filler
    if (!length(.Object@proctime)) .Object@proctime <- proc.time() - tic

    ## return
    .Object
}
