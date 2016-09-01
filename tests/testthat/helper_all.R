data(uciCT)

data <- CharTraj
data_list <- lapply(CharTraj, reinterpolate, newLength = 180)
data_matrix <- do.call(rbind, data_list)
data_subset <- data[1:20]
data_multivariate <- lapply(seq(1L, 100L, 5L), function(x) cbind(data_list[[x]], data_list[[x+1L]]))

labels <- CharTrajLabels
labels_subset <- labels[1:20]
labels_shuffled <- sample(labels)

internal_cvis <- c("Sil", "D", "DB", "DBstar", "CH", "SF", "COP")
external_cvis <- c("RI", "ARI", "J", "FM", "VI")

ctrl <- new("dtwclustControl", window.size = 18L, save.data = FALSE)

## family and proctime will vary from run to run, that's unavoidable
reset_nondeterministic <- function(obj) {
     if (class(obj) == "dtwclust" && validObject(obj)) {
          obj@family <- new("dtwclustFamily")
          obj@proctime[1:5] <- 0
     }

     obj
}

## problems between accuracy of architectures
my_expect_equal_to_reference <- function(object, x32 = FALSE) {
     skip_on_cran()

     if (x32 && grepl("i386", Sys.getenv("R_ARCH")))
          file_name <- paste0("i386/", as.character(substitute(object)), ".rds")
     else
          file_name <- paste0(as.character(substitute(object)), ".rds")

     expect_equal_to_reference(object, file_name)
}
