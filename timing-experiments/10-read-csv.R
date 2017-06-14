path <- "../../CharTrajCSV"
files <- list.files(path)

univariate_labels <- as.factor(sapply(files, substr, start = 1L, stop = 1L))

univariate_series <- lapply(files, function(file) {
    df <- read.csv(file.path(path, file), header = FALSE)
    colnames(df) <- NULL
    dtwclust:::any2list(df)
})

dtwclust:::setnames_inplace(univariate_series, sub(".csv$", "", files))

i <- seq(from = 1L, to = length(univariate_series), by = 3L)
multivariate_series <- lapply(i, function(i) {
    Map(univariate_series[[i]], univariate_series[[i+1L]], univariate_series[[i+2L]],
        f = function(a, b, c) {
            series <- cbind(a, b, c)
            colnames(series) <- sub(".csv$", "", files[i:(i+2L)])
            series
        })
})

dtwclust:::setnames_inplace(multivariate_series, as.character(univariate_labels[i]))

rm("path", "files", "i")
