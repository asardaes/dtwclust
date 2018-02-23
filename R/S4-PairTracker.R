#' Helper for semi-supervised DTW clustering
#'
#' @importFrom methods setRefClass
#'
#' @field xptr External pointer (C++ class). See corresponding file in src/utils/
#'
#' @keywords internal
#'
PairTracker <- methods::setRefClass(
    "PairTracker",
    fields = list(
        xptr = "externalptr"
    ),
    methods = list(
        initialize = function(max_size) {
            "Initialization of C++ helper"
            max_size <- as.integer(max_size)
            if (max_size < 1L) stop("Invalid size")
            # initialize C++ class
            xptr <<- .Call(C_PairTracker__new, max_size, PACKAGE = "dtwclust")
            # return
            invisible(NULL)
        },
        link = function(i, j, link_type) {
            "Link indices i and j.
            Link types: dont_know = -1, cannot_link = 0, must_link = 1.
            Returns TRUE if underlying graph is complete/complete/connected after insertion."
            i <- as.integer(i)[1L]
            j <- as.integer(j)[1L]
            link_type <- as.integer(link_type)[1L]
            if (abs(link_type) > 1L) stop("Invalid link type")
            .Call(C_PairTracker__link, xptr, i, j, link_type, PACKAGE = "dtwclust")
        },
        get_unseen_pair = function() {
            "Get a pair that is not contained in any graph,
            NULL means no unseen pairs left."
            .Call(C_PairTracker__getUnseenPair, xptr, PACKAGE = "dtwclust")
        }
    )
)
