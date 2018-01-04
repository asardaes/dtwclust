# ==================================================================================================
# Miscellaneous
# ==================================================================================================

check_consistency <- function(obj, case, ..., clus_type,
                              diff_lengths = FALSE, cent_missing, cent_char,
                              trace = FALSE, silent = TRUE)
{
    case <- match.arg(case, c("ts", "tslist", "vltslist", "window", "dist", "cent"))

    if (case == "ts") {
        if (!is.numeric(obj)) stop("The series must be numeric")
        if (length(obj) < 1L) stop("The series must have at least one point")
        if (anyNA(obj)) stop("There are missing values in the series")

    } else if (case %in% c("tslist", "vltslist")) {
        if (!is.list(obj)) stop("Oops, data should already be a list by this point...")
        if (length(obj) < 1L) stop("Data is empty")
        if (case == "tslist" && different_lengths(obj)) stop("All series must have the same length")

        sapply(obj, check_consistency, case = "ts", ...)

    } else if (case == "window") {
        if (is.null(obj)) stop("Please provide the 'window.size' parameter")
        if (any(obj <= 0L)) stop("Window width must be larger than 0")

        return(as.integer(obj))

    } else if (case == "dist") {
        if (!is.character(obj) || !pr_DB$entry_exists(obj)) {
            if (silent)
                return(FALSE)
            else
                stop("Please provide the name of a valid distance function registered with the ",
                     "'proxy' package.")
        }

        if (diff_lengths) {
            obj <- tolower(obj)
            if ((obj %in% distances_known) && !(obj %in% distances_difflength))
                stop("Only the following distances are supported for series with different length:\n\t",
                     paste(distances_difflength, collapse = "\t"))
            else if (!(obj %in% distances_known) && trace)
                message("Series have different lengths. ",
                        "Please confirm that the provided distance function supports this.")
        }

        # valid registered distance
        return(TRUE)

    } else if (case == "cent") {
        cent_char <- switch(
            clus_type,
            partitional = {
                if (is.character(obj)) {
                    cent_char <- match.arg(obj, centroids_nonfuzzy)

                    if (diff_lengths &&
                        cent_char %in% centroids_included &&
                        !(cent_char %in% centroids_difflength))
                        stop("Only the following centroids are supported for ",
                             "series with different lengths:\n\t",
                             paste(centroids_difflength, collapse = "\t"))

                } else {
                    force(cent_char)
                }

                # return partitional switch
                cent_char
            },
            fuzzy = {
                if (is.character(obj)) {
                    cent_char <- match.arg(obj, centroids_fuzzy)

                    if (diff_lengths && cent_char == "fcm")
                        stop("Fuzzy c-means does not support series with different length.")

                } else {
                    force(cent_char)
                }

                # return fuzzy switch
                cent_char
            },
            hierarchical =, tadpole = {
                if (is.function(obj))
                    force(cent_char)
                else if (!cent_missing)
                    warning("The 'centroid' argument was provided but it wasn't a function, ",
                            "so it was ignored.",
                            call. = FALSE, immediate. = TRUE)
                else
                    cent_char <- paste("PAM",
                                       switch(clus_type,
                                              hierarchical = "(Hierarchical)",
                                              tadpole = "(TADPole)"))

                # return hierarchical/tadpole switch
                cent_char
            }
        )

        # cent case
        return(cent_char)
    }

    invisible(NULL)
}

# Check if list of series have different length
different_lengths <- function(x) { any(diff(lengths(x)) != 0L) }

# Enlist parameters for do.calls
enlist <- function(..., dots = NULL) { c(list(...), dots) }

# Check if a function has the ellipsis in its formals
has_dots <- function(foo) { is.function(foo) && !is.null(formals(foo)$`...`) }

# Subset dots for do.calls of functions without ellipsis
subset_dots <- function(dots = list(), foo) {
    if (has_dots(foo))
        dots
    else if (length(dots) > 0L)
        dots[intersect(names(dots), names(formals(foo)))]
    else
        list()
}

# Adjust args for tsclust() and TSClusters-class
adjust_args <- function(args, dots) {
    lapply(args, function(arg) {
        arg <- c(arg, dots)
        arg[!duplicated(names(arg))]
    })
}

# Like dynGet() I assume, but that one is supposed to be experimental...
get_from_callers <- function(obj_name, mode = "any") {
    ret <- get0(obj_name, mode = mode, inherits = TRUE)
    if (!is.null(ret)) return(ret)
    for (env in sys.frames()) {
        ret <- get0(obj_name, env, mode = mode, inherits = FALSE)
        if (!is.null(ret)) return(ret)
    }
    stop("Could not find object '", obj_name, "' of mode '", mode, "'") # nocov
}

# ==================================================================================================
# Helper C/C++ functions
# ==================================================================================================

# Create combinations of all possible pairs
call_pairs <- function(n = 2L) {
    n <- as.integer(n)
    if (n < 2L) stop("At least two elements are needed to create pairs.")
    .Call(C_pairs, n, PACKAGE = "dtwclust")
}

# Modify names in place
setnames_inplace <- function(vec, names) {
    if (!is.vector(vec) || !is.vector(names)) stop("Both 'vec' and 'names' must be vectors.")
    if (length(vec) != length(names)) stop("Length mismatch when changing names in place.")
    if (!is.character(names)) stop("Trying to set names in place with non-character names.")
    invisible(.Call(C_setnames_inplace, vec, names, PACKAGE = "dtwclust"))
}

# ==================================================================================================
# Parallel helper functions
# ==================================================================================================

# Custom binary operator for %dopar% to avoid unnecessary warnings
`%op%` <- function(obj, ex) {
    withCallingHandlers({
        ret <- eval.parent(substitute(obj %dopar% ex))
    },
    warning = function(w) {
        if (grepl("package:dtwclust", w$message, ignore.case = TRUE))
            invokeRestart("muffleWarning")
    })

    ret
}

# Split a given object into chunks for parallel workers
#' @importFrom parallel splitIndices
#'
split_parallel <- function(obj, margin = NULL) {
    num_workers <- foreach::getDoParWorkers()
    if (num_workers == 1L) return(structure(list(obj), endpoints = 1L))

    num_tasks <- if (is.null(margin)) length(obj) else dim(obj)[margin]
    if (!is.integer(num_tasks)) stop("Invalid attempt to split an object into parallel tasks")

    num_tasks <- parallel::splitIndices(num_tasks, num_workers)
    num_tasks <- num_tasks[lengths(num_tasks, use.names = FALSE) > 0L]

    if (is.null(margin))
        ret <- lapply(num_tasks, function(id) obj[id])
    else
        ret <- switch(EXPR = margin,
                      lapply(num_tasks, function(id) obj[id, , drop = FALSE]),
                      lapply(num_tasks, function(id) obj[ , id, drop = FALSE]))

    attr(ret, "endpoints") <- lapply(num_tasks, function(ids) { ids[1L] })
    ret
}

# This only works if it's used after split_parallel()
validate_pairwise <- function(x, y) {
    if (!identical(lengths(x, use.names = FALSE), lengths(y, use.names = FALSE)))
        stop("Pairwise distances require the same amount of series in 'x' and 'y'.")

    invisible(NULL)
}

# ==================================================================================================
# Helper distance-related
# ==================================================================================================

# allocate distance matrix for custom proxy loops
#' @importFrom bigmemory big.matrix
#'
allocate_distmat <- function(x_len, y_len, pairwise, symmetric) {
    if (foreach::getDoParWorkers() > 1L) {
        seed <- get0(".Random.seed", .GlobalEnv, mode = "integer") # undo big.matrix() seed change...
        if (pairwise)
            D <- bigmemory::big.matrix(x_len, 1L, "double", 0)
        else if (symmetric)
            D <- bigmemory::big.matrix(x_len, x_len, "double", 0)
        else
            D <- bigmemory::big.matrix(x_len, y_len, "double", 0)
        assign(".Random.seed", seed, .GlobalEnv)

    } else {
        if (pairwise)
            D <- matrix(0, x_len, 1L)
        else if (symmetric)
            D <- matrix(0, x_len, x_len)
        else
            D <- matrix(0, x_len, y_len)
    }
    # return
    D
}

# get endpoints for parallel symmetric distance matrix calculations based on number of workers
symmetric_loop_endpoints <- function(n) {
    if (n < 2L) stop("No symmetric calculations possible for a 1x1 distance matrix")
    num_workers <- foreach::getDoParWorkers()
    if (num_workers == 1L) return(list(list(
        start = list(i = 2L, j = 1L), end = list(i = n, j = n - 1L))
    ))

    # single to double index for symmetric matrices
    s2d <- function(id, n) {
        if (id < n) return(list(i = id + 1L, j = 1L))
        # start at second column
        i <- 3L
        j <- 2L
        start_pair <- n
        end_pair <- n * 2L - 3L
        # j is ready after this while loop finishes
        while (!(id >= start_pair && id <= end_pair)) {
            start_pair <- end_pair + 1L
            end_pair <- start_pair + n - j - 2L
            i <- i + 1L
            j <- j + 1L
        }
        # while loop for i
        while (start_pair < id) {
            i <- i + 1L
            start_pair <- start_pair + 1L
        }
        # return
        list(i = i, j = j)
    }

    num_pairs <- as.integer(n * (n + 1L) / 2L - n)
    if (num_pairs < num_workers) num_workers <- num_pairs
    start_pairs <- cumsum(c(1L, rep(as.integer(num_pairs / num_workers), num_workers - 1L)))
    end_pairs <- c(start_pairs[-1L] - 1L, num_pairs)
    Map(start_pairs, end_pairs, f = function(start_pair, end_pair) {
        list(start = s2d(start_pair, n), end = s2d(end_pair, n))
    })
}

# Euclidean norm
l2norm <- function(x) { sqrt(sum(x * x)) }

# Get the clusinfo slot for TSClusters
compute_clusinfo <- function(k, cluster, cldist) {
    cluster <- factor(cluster, levels = 1L:k)
    # return
    data.frame(
        size = as.integer(table(cluster)),
        av_dist = as.numeric(tapply(cldist[, 1L], list(cluster), mean, default = 0))
    )
}

# ==================================================================================================
# Multivariate helpers
# ==================================================================================================

is_multivariate <- function(x) {
    if (length(x) == 0L) stop("Empty list of series received.")
    ncols <- sapply(x, NCOL)
    if (any(diff(ncols) != 0L)) stop("Inconsistent dimensions across series.")
    any(ncols > 1L)
}

reshape_multivariate <- function(series, cent) {
    ncols <- ncol(series[[1L]])

    series <- lapply(1L:ncols, function(idc) {
        lapply(series, function(s) { s[ , idc, drop = TRUE] })
    })

    cent <- lapply(1L:ncols, function(idc) {
        if (is.null(cent))
            NULL
        else
            cent[ , idc, drop = TRUE]
    })

    list(series = series, cent = cent)
}
