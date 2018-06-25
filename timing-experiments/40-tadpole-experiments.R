# ==================================================================================================
# TADPole experiments
# ==================================================================================================

existing_objects <- ls(all.names = TRUE)

# --------------------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------------------

# Take only one variable of each character
series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
# Length of reinterpolated series
new_length <- median(sapply(series, function(s) { length(s[[1L]]) }))

if (short_experiments) {
    # How many of **each** character to consider
    num_series <- seq(from = 30L, to = 50L, by = 20L)
    # Window sizes to consider
    window_sizes <- seq(from = 40L, to = 80L, by = 20L)
    # Which lower bounds to try
    lbs <- "lbk"
    # Number of evaluations
    times <- 5L

} else {
    # How many of **each** character to consider
    num_series <- seq(from = 10L, to = 50L, by = 10L)
    # Window sizes to consider
    window_sizes <- seq(from = 20L, to = 100L, by = 20L)
    # Which lower bounds to try
    lbs <- c("lbk", "lbi")
    # Number of evaluations
    times <- 30L
}

t1 <- proc.time()

# --------------------------------------------------------------------------------------------------
# Experiments
# --------------------------------------------------------------------------------------------------

# NOTE: all clustering experiments will use tsclust() to include overhead of corresponding checks

mycat("\tRunning TADPole experiments\n")
clus_tadpole_results <- dplyr::bind_rows(lapply(num_series, function(num_series) {
    mycat("\t\t")

    # Get subset and reinterpolate to equal length
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)
    series <- reinterpolate(series, new.length = new_length)

    benchmarks <- dplyr::bind_rows(lapply(lbs, function(lb) {
        benchmark <- lapply(window_sizes, function(window_size) {
            median_time <- median(sapply(1L:times, function(dummy) {
                tsc <- tsclust(series = series, k = 20L, type = "tadpole",
                               trace = FALSE, error.check = FALSE,
                               control = tadpole_control(dc = 10,
                                                         window.size = window_size,
                                                         lb = lb))
                tsc@proctime[["elapsed"]]
            }))

            cat(".")
            # Return data frame
            data.frame(num_series = length(series),
                       k = 20L,
                       dc = 10,
                       window_size = window_size,
                       lb = lb,
                       median_time_s = median_time)
        })

        dplyr::bind_rows(benchmark)
    }))

    cat("\n")
    benchmarks
}))

# Add some metadata
attr(clus_tadpole_results, "proctime") <- proc.time() - t1
attr(clus_tadpole_results, "times") <- times

# ==================================================================================================
# finish
# ==================================================================================================

# Clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "clus_tadpole_results")))
save("clus_tadpole_results", file = "tadpole-results.RData")
cat("\n")
