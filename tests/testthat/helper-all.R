data(uciCT)
options(deparse.max.lines = 2L)

# environment to save objects for regression tests
persistent <- new.env()

# data
data <- CharTraj
data_subset <- data[1L:20L]
data_reinterpolated <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
data_reinterpolated_subset <- data_reinterpolated[1L:20L]
data_multivariate <- CharTrajMV[1L:20L]
data_matrix <- do.call(rbind, data_reinterpolated)

# labels
labels <- CharTrajLabels
labels_subset <- labels[1L:20L]
labels_shuffled <- sample(labels)

# family (environments) and proctime will vary from run to run, that's unavoidable
reset_nondeterministic <- function(obj, clear.data = TRUE) {
    if (inherits(obj, "TSClusters") && methods::validObject(obj)) {
        obj@family <- new("tsclustFamily")
        obj@proctime[1L:5L] <- 0
        if (inherits(obj@control$distmat, "Distmat")) obj@control$distmat <- obj@distmat
        if (clear.data) obj@datalist <- list()
    }
    obj
}

# return file name for rds files
file_name <- function(object, x32 = FALSE) {
    if (x32 && .Machine$sizeof.pointer == 4L)
        paste0("rds/x32/", as.character(substitute(object)), ".rds")
    else
        paste0("rds/", as.character(substitute(object)), ".rds")
}
