data(uciCT)
options(deparse.max.lines = 2L)
suppressWarnings(RNGversion("3.5.0"))
options(dtwclust_sdtw_cent_return_attrs = FALSE,
        dtwclust_suggest_bigmemory = FALSE)

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

# expect_equal for specific slots
expect_equal_slots <- function(current, target, slots = c("cluster", "centroids", "cldist"), ...) {
    for (object_slot in slots) {
        expect_equal(slot(current, object_slot),
                     slot(target, object_slot),
                     ...,
                     info = paste("slot =", object_slot))
    }
}

expect_known_rds <- function(object, path, ..., info = NULL, update = TRUE) {
    file <- if (missing(path)) paste0("rds/", rlang::enexpr(object)) else path

    if (!file.exists(file)) {
        warning("Creating reference value", call. = FALSE)
        saveRDS(object, file, version = 2)
        succeed()
    }
    else {
        ref_val <- readRDS(file)
        comp <- compare(object, ref_val, ...)
        if (update && !comp$equal) {
            saveRDS(object, file, version = version)
        }
        expect(comp$equal,
               sprintf("%s has changed from known value recorded in %s.\n%s",
                       file,
                       encodeString(file, quote = "'"),
                       comp$message),
               info = info)
    }
    invisible(object)
}
