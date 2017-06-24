# ==================================================================================================
# time series prototypes
# ==================================================================================================

existing_objects <- ls(all.names = TRUE)

# --------------------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------------------

#' The dataset has series with different lengths. This threshold specifies the difference there
#' should be between series in the subset that is taken.
length_diff_threshold <- 20L

# Take only one variable of each character
series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
# Get length of each character's time series
len <- sapply(series, function(s) { lengths(s)[1L] })
# Sort lengths
id_ascending <- sort(len, decreasing = FALSE, index.return = TRUE)

# Identify those characters that have length differences greater than the threshold specified
len <- id_ascending$x
id_ascending <- id_ascending$ix
dlen <- diff(len) < length_diff_threshold
while (any(dlen)) {
    rem <- which(dlen)[1L] + 1L
    len <- len[-rem]
    id_ascending <- id_ascending[-rem]
    dlen <- diff(len) < length_diff_threshold
}

# Get the resulting characters
series <- series[id_ascending]
# Also get the multivariate versions
series_mv <- multivariate_series[id_ascending]
# Normalize
series_normalized <- lapply(series, zscore)
series_mv_normalized <- lapply(series_mv, zscore)
# Number of series to consider during prototyping
num_series <- seq(from = 10L, to = 100L, by = 10L)
# Window sizes for the experiments (where applicable)
window_sizes <- seq(from = 10L, to = 100L, by = 10L)
# Number of times each experiment will be repeated (by microbenchmark)
times <- if (short_experiments) 10L else 50L

if (short_experiments) {
    window_sizes <- num_series <- c(10L, 20L)
}

#' NOTE: all experiments are pretty much equivalent. They are run within a new environment so that
#' they don't change variables in the global environment (this one).

# --------------------------------------------------------------------------------------------------
# univariate shape_extraction
# --------------------------------------------------------------------------------------------------

cat("\tRunning shape_extraction experiments for univariate series\n")
cent_shape_univariate <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series_normalized, function(this_series) {
            # Build expressions to evaluate, substituting number of series
            expressions <- lapply(num_series, function(ns) {
                bquote(
                    shape_extraction(this_series[1L:.(ns)], ref, error.check = FALSE)
                )
            })

            # Extract reference series
            ref <- this_series[[length(this_series)]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

            # Return data frame with results
            data.frame(cent = "shape_univariate",
                       series_length = NROW(ref),
                       num_series = num_series,
                       median_time_ms = benchmark$median)
        })

        # Bind results for all series and return to global environment
        plyr::rbind.fill(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# multivariate shape_extraction
# --------------------------------------------------------------------------------------------------

cat("\tRunning shape_extraction experiments for multivariate series\n")
cent_shape_multivariate <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series_mv_normalized, function(this_series) {
            # Build expressions to evaluate, substituting number of series
            expressions <- lapply(num_series, function(ns) {
                bquote(
                    shape_extraction(this_series[1L:.(ns)], ref, error.check = FALSE)
                )
            })

            # Extract reference series
            ref <- this_series[[length(this_series)]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

            # Return data frame with results
            data.frame(cent = "shape_multivariate",
                       series_length = NROW(ref),
                       num_series = num_series,
                       median_time_ms = benchmark$median)
        })

        # Bind results for all series and return to global environment
        plyr::rbind.fill(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# univariate dba
# --------------------------------------------------------------------------------------------------

cat("\tRunning dba experiments for univariate series\n")
cent_dba_univariate <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series, function(this_series) {
            # Build expressions to evaluate, substituting number of series and window size
            expressions <- lapply(num_series, function(ns) {
                lapply(window_sizes, function(ws) {
                    bquote(
                        DBA(this_series[1L:.(ns)], ref, window.size = .(ws), error.check = FALSE)
                    )
                })
            })
            expressions <- unlist(expressions, recursive = FALSE)

            # Extract reference series
            ref <- this_series[[length(this_series)]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

            # Return data frame with results
            data.frame(cent = "dba_univariate",
                       series_length = NROW(ref),
                       window_size = window_sizes,
                       num_series = rep(num_series, each = length(window_sizes)),
                       median_time_ms = benchmark$median)
        })

        # Bind results for all series and return to global environment
        plyr::rbind.fill(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# multivariate dba byS
# --------------------------------------------------------------------------------------------------

cat("\tRunning dba experiments for multivariate series (by-series)\n")
cent_dba_multivariate_byS <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series_mv, function(this_series) {
            # Build expressions to evaluate, substituting number of series and window size
            expressions <- lapply(num_series, function(ns) {
                lapply(window_sizes, function(ws) {
                    bquote(
                        DBA(this_series[1L:.(ns)], ref, window.size = .(ws),
                            mv.ver = "by-s", error.check = FALSE)
                    )
                })
            })
            expressions <- unlist(expressions, recursive = FALSE)

            # Extract reference series
            ref <- this_series[[length(this_series)]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

            # Return data frame with results
            data.frame(cent = "dba_multivariate_byS",
                       series_length = NROW(ref),
                       window_size = window_sizes,
                       num_series = rep(num_series, each = length(window_sizes)),
                       median_time_ms = benchmark$median)
        })

        # Bind results for all series and return to global environment
        plyr::rbind.fill(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# multivariate dba byV
# --------------------------------------------------------------------------------------------------

cat("\tRunning dba experiments for multivariate series (by-variable)\n")
cent_dba_multivariate_byV <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series_mv, function(this_series) {
            # Build expressions to evaluate, substituting number of series and window size
            expressions <- lapply(num_series, function(ns) {
                lapply(window_sizes, function(ws) {
                    bquote(
                        DBA(this_series[1L:.(ns)], ref, window.size = .(ws),
                            mv.ver = "by-v", error.check = FALSE)
                    )
                })
            })
            expressions <- unlist(expressions, recursive = FALSE)

            # Extract reference series
            ref <- this_series[[length(this_series)]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

            # Return data frame with results
            data.frame(cent = "dba_multivariate_byV",
                       series_length = NROW(ref),
                       window_size = window_sizes,
                       num_series = rep(num_series, each = length(window_sizes)),
                       median_time_ms = benchmark$median)
        })

        # Bind results for all series and return to global environment
        plyr::rbind.fill(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# aggregate
# --------------------------------------------------------------------------------------------------

cent_results <- plyr::rbind.fill(
    cent_shape_univariate,
    cent_shape_multivariate,
    cent_dba_univariate,
    cent_dba_multivariate_byS,
    cent_dba_multivariate_byV
)

# Make factor with the given order
cent_results$cent <- factor(cent_results$cent,
                            levels = unique(cent_results$cent))

# Add some metadata
attr(cent_results, "proctime") <- proc.time() - tic
attr(cent_results, "times") <- times

# Clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "cent_results")))

# ==================================================================================================
# finish
# ==================================================================================================

save("cent_results", file = "cent-results.RData")
