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
          if (!is.null(dim(obj)) && min(dim(obj)) != 1) {
               stop("The series must be univariate vectors")
          }
          if (length(obj) < 1) {
               stop("The series must have at least one point")
          }
          if (any(is.na(obj))) {
               stop("There are missing values in the series")
          }

     } else if (case == "tslist") {
          if (class(obj) != "list")
               stop("Series must be provided in the form of a list")

          if (length(obj) < 1)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) {
               ifelse(!is.null(dim(obj)) && min(dim(obj)) != 1, TRUE, FALSE)
          })))
               stop("Each element of the list must be a univariate vector")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (length(unique(sapply(obj, length))) > 1)
               stop("All series must have the same length")

          if (any(sapply(obj, length) < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "vltslist") {

          ## list of variable-length time series

          if (class(obj) != "list")
               stop("Series must be provided in the form of a list")

          if (length(obj) < 1)
               stop("Data is empty")

          if (any(sapply(obj, function(obj) {
               ifelse(!is.null(dim(obj)) && min(dim(obj)) != 1, TRUE, FALSE)
          })))
               stop("Each element of the list must be a univariate vector")

          if (any(!sapply(obj, is.numeric)))
               stop("Each element of the list must be a numeric vector")

          if (any(sapply(obj, length) < 1))
               stop("All series must have at least one element")

          if (any( sapply(obj, function(ts) {any(is.na(ts))}) ))
               stop("Time series cannot have missing elements")

     } else if (case == "window") {
          if (is.null(obj)) {
               stop("Please provide the 'window.size' parameter")
          }
          if (obj <= 0) {
               stop("Window width must be larger than 0")
          }

          return(round(obj))

     } else if (case == "tsmat") {

          if (is.matrix(obj))
               obj <- lapply(seq_len(nrow(obj)), function(i) obj[i,])
          else if (is.numeric(obj))
               obj <- list(obj)
          else if (!is.list(obj))
               stop("Unsupported type for data")

          return(obj)

     } else if (case == "dist") {
          included <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd")
          valid <- c("dtw", "dtw2", "sbd")

          if (is.character(obj) && (obj %in% included) && !(obj %in% valid))
               stop("Only the following distances are supported for series of different lengths:\n\tdtw\tdtw2\tsbd")

     } else if (case == "cent") {
          included <- c("mean", "median", "shape", "dba", "pam")
          valid <- c("dba", "pam", "shape")

          if (is.character(obj) && (obj %in% included) && !(obj %in% valid))
               stop("Only the following centroids are supported for series of different lengths:\n\tdba\tpam\tshape")

     } else {
          stop("Possibly a typo in function consistency_check")
     }

     invisible(NULL)
}

# ========================================================================================================
# Helper functions
# ========================================================================================================

# Create combinations of all possible pairs
call_pairs <- function(n = 2L, byrow = TRUE) {
     if (n < 2)
          stop("I need at least two elements to create pairs.")

     .Call("pairs", n, byrow, PACKAGE = "dtwclust")
}

# Is there a registered parallel backend?
check_parallel <- function(distance = NULL) {
     ret <- foreach::getDoParRegistered() && foreach::getDoParWorkers() > 1L

     if (!is.null(distance))
          ret <- ret && pr_DB$get_entry(distance)$loop

     ret
}

# Split a given number to tasks into parallel workers
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

dtwclust_family <- function(name, dist, allcent) {
     cluster <- function (x, centers, distmat = NULL) {
          if (is.null(distmat))
               distmat <- dist(x, centers)

          max.col(-distmat, "first")
     }

     list(name = name, dist = dist, allcent = allcent, cluster = cluster)
}
