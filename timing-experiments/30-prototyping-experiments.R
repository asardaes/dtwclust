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
window_sizes <- seq(from = 30L, to = 90L, by = 30L)
# Number of times each experiment will be repeated (by microbenchmark)
times <- 30L
# Number of threads to test for those that support multi-threading
num_threads <- c(1L, 2L, 4L)

if (short_experiments) {
    times <- 10L
    window_sizes <- num_series <- c(10L, 20L)
    num_threads <- c(2L, 4L)
}

#' NOTE: all experiments are pretty much equivalent. They are run within a new environment so that
#' they don't change variables in the global environment (this one).

# For metadata
t1 <- proc.time()

# --------------------------------------------------------------------------------------------------
# univariate shape_extraction
# --------------------------------------------------------------------------------------------------

mycat("\tRunning shape_extraction experiments for univariate series\n")
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
        dplyr::bind_rows(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# multivariate shape_extraction
# --------------------------------------------------------------------------------------------------

mycat("\tRunning shape_extraction experiments for multivariate series\n")
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
        dplyr::bind_rows(benchmarks)
    })

# --------------------------------------------------------------------------------------------------
# univariate dba
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dba experiments for univariate series\n")
cent_dba_univariate <- with(
    new.env(),
    {
        # Loop across number of threads to test
        dplyr::bind_rows(lapply(num_threads, function(num_threads) {
            mycat("\t\t")
            RcppParallel::setThreadOptions(num_threads)

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
                cat(".")
                data.frame(cent = "dba_univariate",
                           num_threads = num_threads,
                           series_length = NROW(ref),
                           window_size = window_sizes,
                           num_series = rep(num_series, each = length(window_sizes)),
                           median_time_ms = benchmark$median)
            })

            cat("\n")
            # Bind results for all series
            dplyr::bind_rows(benchmarks)
        }))
    })

# --------------------------------------------------------------------------------------------------
# multivariate dba byS
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dba experiments for multivariate series (by-series)\n")
cent_dba_multivariate_byS <- with(
    new.env(),
    {
        # Loop across number of threads to test
        dplyr::bind_rows(lapply(num_threads, function(num_threads) {
            mycat("\t\t")
            RcppParallel::setThreadOptions(num_threads)

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
                cat(".")
                data.frame(cent = "dba_multivariate_byS",
                           num_threads = num_threads,
                           series_length = NROW(ref),
                           window_size = window_sizes,
                           num_series = rep(num_series, each = length(window_sizes)),
                           median_time_ms = benchmark$median)
            })

            cat("\n")
            # Bind results for all series and return to global environment
            dplyr::bind_rows(benchmarks)
        }))
    })

# --------------------------------------------------------------------------------------------------
# multivariate dba byV
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dba experiments for multivariate series (by-variable)\n")
cent_dba_multivariate_byV <- with(
    new.env(),
    {
        # Loop across number of threads to test
        dplyr::bind_rows(lapply(num_threads, function(num_threads) {
            mycat("\t\t")
            RcppParallel::setThreadOptions(num_threads)

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
                cat(".")
                data.frame(cent = "dba_multivariate_byV",
                           num_threads = num_threads,
                           series_length = NROW(ref),
                           window_size = window_sizes,
                           num_series = rep(num_series, each = length(window_sizes)),
                           median_time_ms = benchmark$median)
            })

            cat("\n")
            # Bind results for all series and return to global environment
            dplyr::bind_rows(benchmarks)
        }))
    })

# --------------------------------------------------------------------------------------------------
# univariate sdtw_cent
# --------------------------------------------------------------------------------------------------

mycat("\tRunning sdtw_cent experiments for univariate series\n")
cent_sdtw_univariate <- with(
    new.env(),
    {
        # Loop across number of threads to test
        dplyr::bind_rows(lapply(num_threads, function(num_threads) {
            mycat("\t\t")
            RcppParallel::setThreadOptions(num_threads)

            # Loop along subset of extracted subsets
            series <- series[seq(from = 1L, to = length(series), by = 2L)]
            benchmarks <- lapply(series, function(this_series) {
                # Build expressions to evaluate, substituting number of series
                expressions <- lapply(num_series, function(ns) {
                    bquote(
                        sdtw_cent(this_series[1L:.(ns)], ref, error.check = FALSE, num_threads = num_threads)
                    )
                })

                # Extract reference series
                ref <- this_series[[length(this_series)]]

                # Evaluate expressions
                benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

                # Return data frame with results
                cat(".")
                data.frame(cent = "sdtw_cent_univariate",
                           num_threads = num_threads,
                           series_length = NROW(ref),
                           num_series = num_series,
                           median_time_ms = benchmark$median)
            })

            cat("\n")
            # Bind results for all series and return to global environment
            dplyr::bind_rows(benchmarks)
        }))
    })

# --------------------------------------------------------------------------------------------------
# multivariate sdtw_cent
# --------------------------------------------------------------------------------------------------

mycat("\tRunning sdtw_cent experiments for multivariate series\n")
cent_sdtw_multivariate <- with(
    new.env(),
    {
        # Loop across number of threads to test
        dplyr::bind_rows(lapply(num_threads, function(num_threads) {
            mycat("\t\t")
            RcppParallel::setThreadOptions(num_threads)

            # Loop along subset of extracted subsets
            series_mv <- series_mv[seq(from = 1L, to = length(series_mv), by = 2L)]
            benchmarks <- lapply(series_mv, function(this_series) {
                # Build expressions to evaluate, substituting number of series
                expressions <- lapply(num_series, function(ns) {
                    bquote(
                        sdtw_cent(this_series[1L:.(ns)], ref, error.check = FALSE, num_threads = num_threads)
                    )
                })

                # Extract reference series
                ref <- this_series[[length(this_series)]]

                # Evaluate expressions
                benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

                # Return data frame with results
                cat(".")
                data.frame(cent = "sdtw_cent_multivariate",
                           num_threads = num_threads,
                           series_length = NROW(ref),
                           num_series = num_series,
                           median_time_ms = benchmark$median)
            })

            cat("\n")
            # Bind results for all series and return to global environment
            dplyr::bind_rows(benchmarks)
        }))
    })

# --------------------------------------------------------------------------------------------------
# aggregate
# --------------------------------------------------------------------------------------------------

cent_results <- dplyr::bind_rows(
    cent_shape_univariate,
    cent_shape_multivariate,
    cent_dba_univariate,
    cent_dba_multivariate_byS,
    cent_dba_multivariate_byV,
    cent_sdtw_univariate,
    cent_sdtw_multivariate
)

# Make factor with the given order
cent_results$cent <- factor(cent_results$cent,
                            levels = unique(cent_results$cent))

# Add some metadata
attr(cent_results, "proctime") <- proc.time() - t1
attr(cent_results, "times") <- times

# ==================================================================================================
# finish
# ==================================================================================================

# Clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "cent_results")))
save("cent_results", file = "cent-results.RData")
cat("\n")
