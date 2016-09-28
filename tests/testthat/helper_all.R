data(uciCT)

data <- CharTraj
data_subset <- data[1L:20L]
data_reinterpolated <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
data_multivariate <- lapply(seq(1L, 100L, 5L), function(x) {
     cbind(data_reinterpolated[[x]], data_reinterpolated[[x+1L]])
})
data_matrix <- do.call(rbind, data_reinterpolated)

labels <- CharTrajLabels
labels_subset <- labels[1L:20L]
labels_shuffled <- sample(labels)

internal_cvis <- c("Sil", "D", "DB", "DBstar", "CH", "SF", "COP")
external_cvis <- c("RI", "ARI", "J", "FM", "VI")

## family (environments) and proctime will vary from run to run, that's unavoidable
reset_nondeterministic <- function(obj, clear.data = TRUE) {
     if (class(obj) == "dtwclust" && validObject(obj)) {
          obj@family <- new("dtwclustFamily")
          obj@proctime[1:5] <- 0

          if (clear.data) obj@datalist <- list()
     }

     obj
}

file_name <- function(object) {
     paste0("rds/", as.character(substitute(object)), ".rds")
}
