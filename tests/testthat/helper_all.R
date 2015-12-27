data(uciCT)

data <- CharTraj
data_list <- lapply(CharTraj, reinterpolate, newLength = 180)
data_matrix <- do.call(rbind, data_list)
data_subset <- data[1:20]

labels <- CharTrajLabels
labels_subset <- labels[1:20]

ctrl <- new("dtwclustControl", window.size = 18L)

## family and proctime will vary from run to run, that's unavoidable
reset_nondeterministic <- function(obj) {
     if (class(obj) == "dtwclust" && validObject(obj)) {
          obj@family <- new("dtwclustFamily")
          obj@proctime[1:5] <- 0
     }

     obj
}
