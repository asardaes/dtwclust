# ========================================================================================================
# Miscellaneous
# ========================================================================================================

check_consistency <- function(obj, case, ..., trace = FALSE, Lengths = FALSE, silent = TRUE) {
    case <- match.arg(case, c("ts", "tslist", "vltslist", "window", "dist", "cent"))

    if (case == "ts") {
        if (!is.numeric(obj)) {
            stop("The series must be numeric")
        }
        if (length(obj) < 1L) {
            stop("The series must have at least one point")
        }
        if (any(is.na(obj))) {
            stop("There are missing values in the series")
        }

    } else if (case %in% c("tslist", "vltslist")) {
        if (!is.list(obj))
            stop("Oops, data should already be a list by this point...")

        if (length(obj) < 1L)
            stop("Data is empty")

        if (case == "tslist" && different_lengths(obj))
            stop("All series must have the same length")

        sapply(obj, check_consistency, case = "ts", ...)

    } else if (case == "window") {
        if (is.null(obj))
            stop("Please provide the 'window.size' parameter")

        if (obj <= 0L)
            stop("Window width must be larger than 0")

        return(as.integer(obj))

    } else if (case == "dist") {
        included <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd", "dtw_basic", "gak")
        valid <- c("dtw", "dtw2", "sbd", "dtw_basic", "gak") ## for different lengths

        if (!is.character(obj) || !pr_DB$entry_exists(obj)) {
            if (silent)
                return(FALSE)
            else
                stop("Please provide the name of a valid distance function registered with the ",
                     "'proxy' package.")
        }

        if (Lengths) {
            if ((obj %in% included) && !(obj %in% valid))
                stop("Only the following distances are supported for series with different length:\n\t",
                     paste(valid, collapse = "\t"))
            else if (!(obj %in% included) && trace)
                message("Series have different lengths. Please confirm that the provided distance ",
                        "function supports this.")
        }

        ## valid registered distance
        return(TRUE)

    } else if (case == "cent") {
        ## only checking for different lengths
        included <- c("mean", "median", "shape", "dba", "pam", "fcm", "fcmdd")
        valid <- c("dba", "pam", "shape", "fcmdd")

        if (is.character(obj) && (obj %in% included) && !(obj %in% valid))
            stop("Only the following centroids are supported for series with different lengths:",
                 "\n\tdba\tpam\tshape\tfcmdd")

    }

    invisible(NULL)
}

# Coerce to list
any2list <- function(obj) {
    if (is.matrix(obj)) {
        Names <- rownames(obj)
        obj <- lapply(seq_len(nrow(obj)), function(i) obj[i, ])
        names(obj) <- Names

    } else if (is.numeric(obj)) {
        obj <- list(obj)

    } else if (is.data.frame(obj)) {
        obj <- any2list(as.matrix(obj))

    } else if (!is.list(obj))
        stop("Unsupported data type.")

    obj
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

# This only works if it's used after split_parallel()
validate_pairwise <- function(x, y) {
    if (!identical(lengths(x, use.names = FALSE), lengths(y, use.names = FALSE)))
        stop("Pairwise distances require the same amount of series in 'x' and 'y'.")

    invisible(NULL)
}

# ========================================================================================================
# Helper C/C++ functions
# ========================================================================================================

# Create combinations of all possible pairs
call_pairs <- function(n = 2L, lower = TRUE) {
    if (n < 2L) stop("At least two elements are needed to create pairs.")

    .Call("pairs", n, lower, PACKAGE = "dtwclust")
}

# ========================================================================================================
# Parallel helper functions
# ========================================================================================================

# allow_parallel <- function() {
#     do_par <- TRUE
#
#     function(flag) {
#         if (!missing(flag))
#             do_par <<- as.logical(flag)[1L]
#
#         invisible(do_par)
#     }
# }
#
# Allow parallel computation with \pkg{dtwclust} functions
#
# Set an internal flag that is checked by the functions that support parallel computation and
# determines if they should attempt to use a parallel backend if one is registered.
#
# @export
#
# @param flag \code{TRUE} to allow use of parallel backends or \code{FALSE} to prevent it.
#
# @return \code{flag} invisibly
#
# parallel_dtwclust <- allow_parallel()

# Custom binary operator for %dopar% to avoid unnecessary warnings
`%op%` <- function(obj, ex) {
    withCallingHandlers({
        ret <- eval.parent(substitute(obj %dopar% ex))

    }, warning = function(w) {
        if (grepl("package:dtwclust", w$message, ignore.case = TRUE))
            invokeRestart("muffleWarning")
    })

    ret
}

# Split a given object into chunks for parallel workers
split_parallel <- function(obj, margin = NULL) {
    num_workers <- foreach::getDoParWorkers()

    if (num_workers == 1L)
        return(list(obj))

    if (is.null(margin))
        num_tasks <- length(obj)
    else
        num_tasks <- dim(obj)[margin]

    if (is.na(num_tasks))
        stop("Invalid attempt to split an object into parallel tasks")

    num_tasks <- parallel::splitIndices(num_tasks, num_workers)
    num_tasks <- num_tasks[lengths(num_tasks, use.names = FALSE) > 0L]

    if (is.null(margin))
        ret <- lapply(num_tasks, function(id) obj[id])
    else
        ret <- switch(EXPR = margin,
                      lapply(num_tasks, function(id) obj[id, , drop = FALSE]),
                      lapply(num_tasks, function(id) obj[ , id, drop = FALSE])
        )

    ret
}

## tasks created based on getDoParWorkers() could be larger than tasks based on objects
allocate_matrices <- function(mat = NULL, ..., target.size) {
    if (foreach::getDoParWorkers() > 1L) {
        MAT <- lapply(1L:target.size, function(dummy) {
            matrix(0, ...)
        })

    } else if (is.null(mat)) {
        MAT <- list(matrix(0, ...))

    } else {
        MAT <- list(mat)
    }

    MAT
}

# ========================================================================================================
# Helper distance-related
# ========================================================================================================

# column-wise medians
colMedians <- function(mat) { apply(mat, 2L, stats::median) }

# l-norm, I think I'll only use L2
lnorm <- function(x, n = 2) {
    if (n %% 2 == 0L)
        sum(x ^ n) ^ (1 / n)
    else
        sum(abs(x) ^ n) ^ (1 / n)
}

# PREFUN for some of my proxy distances so that they support 'pairwise' directly
proxy_prefun <- function(x, y, pairwise, params, reg_entry) {
    params$pairwise <- pairwise

    list(x = x, y = y,
         pairwise = pairwise,
         p = params,
         reg_entry = reg_entry)
}

# ========================================================================================================
# Multviariate helpers
# ========================================================================================================

is_multivariate <- function(x) {
    ncols <- sapply(x, NCOL)

    if (any(diff(ncols) != 0L))
        stop("Inconsistent dimensions across series.")

    any(ncols > 1L)
}

reshape_multviariate <- function(series, cent) {
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
