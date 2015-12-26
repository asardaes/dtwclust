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
          if (length(obj) < 1) {
               stop("The series must have at least one point")
          }
          if (any(is.na(obj))) {
               stop("There are missing values in the series")
          }

     } else if (case == "tslist") {
          Lengths <- sapply(obj, length)

          if (class(obj) != "list")
               stop("Series must be provided in a list")

          if (length(obj) < 1)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) !is.null(dim(obj)) )))
               stop("Each element of the list must be a univariate vector (at least one has non-NULL dim)")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (length(unique(Lengths)) > 1)
               stop("All series must have the same length")

          if (any(Lengths < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "vltslist") {
          Lengths <- sapply(obj, length)

          ## list of variable-length time series

          if (class(obj) != "list")
               stop("Series must be provided in a list")

          if (length(obj) < 1)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) !is.null(dim(obj)) )))
               stop("Each element of the list must be a univariate vector (at least one has non-NULL dim)")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (any(Lengths < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "window") {
          if (is.null(obj))
               stop("Please provide the 'window.size' parameter")

          if (obj <= 0)
               stop("Window width must be larger than 0")

          return(as.integer(obj))

     } else if (case == "tsmat") {

          if (is.matrix(obj)) {
               Names <- rownames(obj)
               obj <- lapply(seq_len(nrow(obj)), function(i) obj[i,])
               names(obj) <- Names

          } else if (is.numeric(obj)) {
               obj <- list(obj)

          } else if (!is.list(obj))
               stop("Unsupported type for data")

          return(obj)

     } else if (case == "dist") {
          .local <- function(obj, trace = FALSE, lengths = FALSE, silent = FALSE, ...) {
               included <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd")
               valid <- c("dtw", "dtw2", "sbd")

               ## pr_DB$entry_exists is acting weird
               if (!is.character(obj) || class(try(pr_DB[[obj]], TRUE)) == "try-error") {
                    if (silent)
                         return(FALSE)
                    else
                         stop("Please provide a valid distance function registered with the proxy package.")
               }

               if ((obj %in% included) && !(obj %in% valid) && lengths)
                    stop("Only the following distances are supported for series with different lengths:\n\tdtw\tdtw2\tsbd")
               else if(!(obj %in% included) && trace && lengths)
                    message("Series have different lengths. Please confirm that the provided distance function supports this.")

               TRUE # valid distance
          }

          return(.local(obj, ...))

     } else if (case == "cent") {
          included <- c("mean", "median", "shape", "dba", "pam")
          valid <- c("dba", "pam", "shape")

          if (is.character(obj) && (obj %in% included) && !(obj %in% valid))
               stop("Only the following centroids are supported for series with different lengths:\n\tdba\tpam\tshape")
          else if(is.character(obj) && !(obj %in% included) && list(...)$trace)
               message("Series have different lengths. Please confirm that the provided centroid function supports this.")

     } else {
          stop("Possibly a typo in function consistency_check")
     }

     invisible(NULL)
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
call_pairs <- function(n = 2L, byrow = TRUE) {
     if (n < 2)
          stop("At least two elements are needed to create pairs.")

     .Call("pairs", n, byrow, PACKAGE = "dtwclust")
}

# Is there a registered parallel backend?
check_parallel <- function(distance = NULL) {
     ret <- foreach::getDoParRegistered() && foreach::getDoParWorkers() > 1L

     if (!is.null(distance))
          ret <- ret && pr_DB$get_entry(distance)$loop

     ret
}

# Split a given object into tasks for parallel workers
split_parallel <- function(obj, tasks, margin = NULL) {
     tasks <- parallel::splitIndices(tasks, foreach::getDoParWorkers())
     tasks <- tasks[sapply(tasks, length, USE.NAMES = FALSE) > 0]

     if (is.null(margin))
          ret <- lapply(tasks, function(id) obj[id])
     else
          ret <- switch(EXPR = margin,
                        lapply(tasks, function(id) obj[id,]),
                        lapply(tasks, function(id) obj[,id])
          )

     ret
}

# column-wise medians
colMedians <- function(mat) { apply(mat, 2, stats::median) }
