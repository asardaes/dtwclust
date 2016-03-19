# ========================================================================================================
# Check consistency, used by other functions
# ========================================================================================================

consistency_check <- function(obj, case, ...) {

     case <- match.arg(case, c("ts", "tslist", "vltslist",
                               "window", "tsmat",
                               "dist", "cent"))

     if (case == "ts") {
          if (!is.numeric(obj)) {
               stop("The series must be numeric")
          }
          if (!is.null(dim(obj))) {
               stop("The series must be univariate vectors (provided object has non-NULL dim)")
          }
          if (length(obj) < 1L) {
               stop("The series must have at least one point")
          }
          if (any(is.na(obj))) {
               stop("There are missing values in the series")
          }

     } else if (case == "tslist") {
          Lengths <- lengths(obj)

          if (class(obj) != "list")
               stop("Series must be provided in a list")

          if (length(obj) < 1L)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) !is.null(dim(obj)) )))
               stop("Each element of the list must be a univariate vector (at least one has non-NULL dim)")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (length(unique(Lengths)) > 1L)
               stop("All series must have the same length")

          if (any(Lengths < 1L))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "vltslist") {
          Lengths <- lengths(obj)

          ## list of variable-length time series

          if (class(obj) != "list")
               stop("Series must be provided in a list")

          if (length(obj) < 1L)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) !is.null(dim(obj)) )))
               stop("Each element of the list must be a univariate vector (at least one has non-NULL dim)")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (any(Lengths < 1L))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "window") {
          if (is.null(obj))
               stop("Please provide the 'window.size' parameter")

          if (obj <= 0L)
               stop("Window width must be larger than 0")

          return(as.integer(obj))

     } else if (case == "tsmat") {

          if (is.matrix(obj)) {
               Names <- rownames(obj)
               obj <- lapply(seq_len(nrow(obj)), function(i) obj[i, ])
               names(obj) <- Names

          } else if (is.numeric(obj)) {
               obj <- list(obj)

          } else if (is.data.frame(obj)) {
               obj <- as.list(obj)

          } else if (!is.list(obj))
               stop("Unsupported type for data")

          return(obj)

     } else if (case == "dist") {
          .local <- function(obj, trace = FALSE, Lengths = FALSE, silent = TRUE, ...) {
               included <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd")
               valid <- c("dtw", "dtw2", "sbd")

               ## pr_DB$entry_exists is acting weird
               if (!is.character(obj) || class(try(pr_DB[[obj]], TRUE)) == "try-error") {
                    if (silent)
                         return(FALSE)
                    else
                         stop("Please provide a valid distance function registered with the proxy package.")
               }

               if (Lengths) {
                    if ((obj %in% included) && !(obj %in% valid))
                         stop("Only the following distances are supported for series with different length:\n\tdtw\tdtw2\tsbd")
                    else if(!(obj %in% included) && trace)
                         message("Series have different length. Please confirm that the provided distance function supports this.")
               }

               TRUE # valid registered distance
          }

          return(.local(obj, ...))

     } else if (case == "cent") {
          included <- c("mean", "median", "shape", "dba", "pam", "fcm")
          valid <- c("dba", "pam", "shape")

          if (is.character(obj) && (obj %in% included) && !(obj %in% valid))
               stop("Only the following centroids are supported for series with different length:\n\tdba\tpam\tshape")
          else if(is.character(obj) && !(obj %in% included) && list(...)$trace)
               message("Series have different length. Please confirm that the provided centroid function supports this.")

     } else {
          stop("Possibly a typo in function consistency_check")
     }

     invisible(NULL)
}

# ========================================================================================================
# Membership update for fuzzy c-means clustering
# ========================================================================================================

fcm_cluster <- function(distmat, m) {
     cprime <- apply(distmat, 1L, function(dist_row) { sum( (1 / dist_row) ^ (2 / (m - 1)) ) })

     u <- 1 / apply(distmat, 2L, function(dist_col) { cprime * dist_col ^ (2 / (m - 1)) })

     if (is.null(dim(u))) u <- rbind(u) # for predict generic

     u
}

# ========================================================================================================
# Fuzzy objective function
# ========================================================================================================

fuzzy_objective <- function(u, distmat, m) {
     sum(u^m * distmat^2)
}

# ========================================================================================================
# Helper functions
# ========================================================================================================

# Running extremes
call_runminmax <- function(series, window) {
     series <- as.numeric(series)

     consistency_check(series, "ts")
     window <- consistency_check(window, "window")

     .Call("runminmax", series, window, PACKAGE = "dtwclust")
}

# Create combinations of all possible pairs
call_pairs <- function(n = 2L, lower = TRUE) {
     if (n < 2L)
          stop("At least two elements are needed to create pairs.")

     .Call("pairs", n, lower, PACKAGE = "dtwclust")
}

# Is there a registered parallel backend?
check_parallel <- function() {
     if (is.null(foreach::getDoParName()))
          foreach::registerDoSEQ()

     foreach::getDoParWorkers() > 1L
}

# Split a given object into tasks for parallel workers
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

# column-wise medians
colMedians <- function(mat) { apply(mat, 2L, stats::median) }

# PREFUN for some of my proxy distances so that they support 'pairwise' direclty
# Currently results in warnings because proxy does !is.na(reg_entry$PREFUN)
proxy_prefun <- function(x, y, pairwise, params, reg_entry) {
     params$pairwise <- pairwise

     list(x = x, y = y,
          pairwise = pairwise,
          p = params,
          reg_entry = reg_entry)
}
