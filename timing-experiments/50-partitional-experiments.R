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
    num_series <- seq(from = 10L, to = 30L, by = 10L)
    # Window sizes to consider
    window_sizes <- seq(from = 40L, to = 60L, by = 20L)
    # Number of evaluations
    times <- 2L

} else {
    # How many of **each** character to consider
    num_series <- seq(from = 5L, to = 30L, by = 5L)
    # Window sizes to consider
    window_sizes <- seq(from = 20L, to = 100L, by = 20L)
    # Number of evaluations
    times <- 30L
}

t1 <- proc.time()

# --------------------------------------------------------------------------------------------------
# dtw_basic vs dtw_lb (PAM)
# --------------------------------------------------------------------------------------------------

# NOTE: all clustering experiments will use tsclust() to include overhead of corresponding checks

cat("\tRunning dtw_basic vs dtw_lb clustering experiments (PAM)\n")
clus_dtwb_dtwlb_pam_results <- plyr::rbind.fill(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset and reinterpolate to equal length
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)
    series <- reinterpolate(series, new.length = new_length)

    benchmarks <- plyr::rbind.fill(lapply(window_sizes, function(window_size) {
        times <- sapply(1L:times, function(dummy) {
            tsc_dtwb <- tsclust(series = series, k = 20L, type = "partitional",
                                distance = "dtw_basic", centroid = "pam",
                                seed = window_size, trace = FALSE, error.check = FALSE,
                                control = partitional_control(pam.precompute = TRUE,
                                                              iter.max = 10L),
                                args = tsclust_args(dist = list(window.size = window_size,
                                                                norm = "L1",
                                                                step.pattern = symmetric1)))

            tsc_dtwbs <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_basic", centroid = "pam",
                                 seed = window_size, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(pam.precompute = FALSE,
                                                               pam.sparse = TRUE,
                                                               iter.max = 10L),
                                 args = tsclust_args(dist = list(window.size = window_size,
                                                                 norm = "L1",
                                                                 step.pattern = symmetric1)))

            tsc_dtwlb <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_lb", centroid = "pam",
                                 seed = window_size, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(pam.precompute = FALSE,
                                                               iter.max = 10L),
                                 args = tsclust_args(dist = list(window.size = window_size,
                                                                 norm = "L1",
                                                                 step.pattern = symmetric1)))

            c(dtwb = tsc_dtwb@proctime[["elapsed"]],
              dtwbs = tsc_dtwbs@proctime[["elapsed"]],
              dtwlb = tsc_dtwlb@proctime[["elapsed"]])
        })

        median_times <- apply(times, 1L, median)

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = 20L,
                   window_size = window_size,
                   dtw_basic_median_time_s = median_times[["dtwb"]],
                   dtw_basic_sparse_median_time_s = median_times[["dtwbs"]],
                   dtw_lb_median_time_s = median_times[["dtwlb"]])
    }))

    cat("\n")
    benchmarks
}))

# Add some metadata
attr(clus_dtwb_dtwlb_pam_results, "proctime") <- proc.time() - t1
attr(clus_dtwb_dtwlb_pam_results, "times") <- times

# --------------------------------------------------------------------------------------------------
# dtw_basic vs dtw_lb (DBA)
# --------------------------------------------------------------------------------------------------

# NOTE: all clustering experiments will use tsclust() to include overhead of corresponding checks

t1 <- proc.time()
cat("\tRunning dtw_basic vs dtw_lb clustering experiments (DBA)\n")
clus_dtwb_dtwlb_dba_results <- plyr::rbind.fill(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset and reinterpolate to equal length
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)
    series <- reinterpolate(series, new.length = new_length)

    benchmarks <- plyr::rbind.fill(lapply(window_sizes, function(window_size) {
        times <- sapply(1L:times, function(dummy) {
            tsc_dtwb <- tsclust(series = series, k = 20L, type = "partitional",
                                distance = "dtw_basic", centroid = "dba",
                                seed = window_size, trace = FALSE, error.check = FALSE,
                                control = partitional_control(iter.max = 10L),
                                args = tsclust_args(dist = list(window.size = window_size,
                                                                norm = "L1",
                                                                step.pattern = symmetric1),
                                                    cent = list(window.size = window_size,
                                                                max.iter = 15L,
                                                                step.pattern = symmetric1)))

            tsc_dtwlb <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_lb", centroid = "dba",
                                 seed = window_size, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(iter.max = 10L),
                                 args = tsclust_args(dist = list(window.size = window_size,
                                                                 norm = "L1",
                                                                 step.pattern = symmetric1),
                                                     cent = list(window.size = window_size,
                                                                 max.iter = 15L,
                                                                 step.pattern = symmetric1)))

            c(dtwb = tsc_dtwb@proctime[["elapsed"]],
              dtwlb = tsc_dtwlb@proctime[["elapsed"]])
        })

        median_times <- apply(times, 1L, median)

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = 20L,
                   window_size = window_size,
                   dtw_basic_median_time_s = median_times[["dtwb"]],
                   dtw_lb_median_time_s = median_times[["dtwlb"]])
    }))

    cat("\n")
    benchmarks
}))

# Add some metadata
attr(clus_dtwb_dtwlb_dba_results, "proctime") <- proc.time() - t1
attr(clus_dtwb_dtwlb_dba_results, "times") <- times

# ==================================================================================================
# finish
# ==================================================================================================

# Clean
rm(list = setdiff(ls(all.names = TRUE),
                  c(existing_objects,
                    "clus_dtwb_dtwlb_pam_results",
                    "clus_dtwb_dtwlb_dba_results")))
save("clus_dtwb_dtwlb_pam_results",
     "clus_dtwb_dtwlb_dba_results",
     file = "partitional-results.RData")
cat("\n")
